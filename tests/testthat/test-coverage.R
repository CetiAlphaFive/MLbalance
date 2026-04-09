# test-coverage.R
# Additional tests to improve code coverage across all source files.
# Targets untested branches, internal helpers, and S3 method edge cases.

# ============================================================================
# UTILS: .make_ci
# ============================================================================

test_that(".make_ci returns correct confidence interval", {
  ci <- MLbalance:::.make_ci(1.0, 0.5, alpha = 0.05)
  expect_length(ci, 2)
  expect_true(ci[1] < 1.0)
  expect_true(ci[2] > 1.0)
  # symmetric around estimate

  expect_equal(1.0 - ci[1], ci[2] - 1.0, tolerance = 1e-10)
})

test_that(".make_ci works with custom alpha", {
  ci_90 <- MLbalance:::.make_ci(0, 1, alpha = 0.10)
  ci_99 <- MLbalance:::.make_ci(0, 1, alpha = 0.01)
  # 99% CI should be wider than 90% CI
  expect_true(ci_99[2] - ci_99[1] > ci_90[2] - ci_90[1])
})

# ============================================================================
# UTILS: .g_theme
# ============================================================================

test_that(".g_theme returns a ggplot theme", {
  th <- MLbalance:::.g_theme()
  expect_s3_class(th, "theme")
})

test_that(".g_theme applies subtitle_color when provided", {
  th <- MLbalance:::.g_theme(subtitle_color = "red")
  expect_s3_class(th, "theme")
  # subtitle element should be present
  expect_true(!is.null(th$plot.subtitle))
})

# ============================================================================
# UTILS: .save_rng_state / .restore_rng_state
# ============================================================================

test_that("RNG state save/restore round-trips", {
  set.seed(42)
  saved <- MLbalance:::.save_rng_state()
  x1 <- runif(1)
  MLbalance:::.restore_rng_state(saved)
  x2 <- runif(1)
  expect_equal(x1, x2)
})

test_that(".restore_rng_state(NULL) removes .Random.seed", {
  set.seed(1)  # ensure .Random.seed exists
  MLbalance:::.restore_rng_state(NULL)
  expect_false(exists(".Random.seed", envir = globalenv()))
  # Re-seed so subsequent tests aren't affected
  set.seed(1995)
})

# ============================================================================
# UTILS: .prepare_covariates_for_grf
# ============================================================================

test_that(".prepare_covariates_for_grf handles numeric-only data frame", {
  df <- data.frame(a = 1:5, b = rnorm(5))
  result <- MLbalance:::.prepare_covariates_for_grf(df)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
})

test_that(".prepare_covariates_for_grf one-hot encodes unordered factors", {
  df <- data.frame(a = 1:6, b = factor(c("x", "y", "z", "x", "y", "z")))
  result <- MLbalance:::.prepare_covariates_for_grf(df)
  expect_true(is.matrix(result))
  # 1 numeric col + 3 dummy cols (model.matrix ~ . - 1 for 3 levels)
  expect_true(ncol(result) >= 4)
})

test_that(".prepare_covariates_for_grf converts ordered factors to numeric", {
  df <- data.frame(
    a = 1:5,
    b = ordered(c("low", "med", "high", "low", "med"), levels = c("low", "med", "high"))
  )
  result <- MLbalance:::.prepare_covariates_for_grf(df)
  expect_true(is.matrix(result))
  # Ordered factor → numeric, so no expansion
  expect_equal(ncol(result), 2)
})

test_that(".prepare_covariates_for_grf converts character columns to factor first", {
  df <- data.frame(a = 1:4, b = c("cat", "dog", "cat", "dog"),
                   stringsAsFactors = FALSE)
  result <- MLbalance:::.prepare_covariates_for_grf(df)
  expect_true(is.matrix(result))
  expect_true(ncol(result) >= 3)  # a + cat + dog dummies
})

# ============================================================================
# UTILS: .permute_treatment
# ============================================================================

test_that(".permute_treatment shuffles without clusters or blocks", {
  set.seed(1)
  T_orig <- rep(0:1, each = 50)
  T_perm <- MLbalance:::.permute_treatment(T_orig)
  expect_equal(length(T_perm), length(T_orig))
  expect_equal(sort(T_perm), sort(T_orig))
  # Very unlikely to be identical
  expect_false(all(T_perm == T_orig))
})

test_that(".permute_treatment respects clusters", {
  set.seed(1)
  clusters <- rep(1:10, each = 10)
  T_orig <- rep(rep(0:1, each = 5), each = 10)
  T_perm <- MLbalance:::.permute_treatment(T_orig, clusters = clusters)
  # Treatment should be constant within clusters after permutation
  for (cl in unique(clusters)) {
    vals <- T_perm[clusters == cl]
    expect_true(length(unique(vals)) == 1)
  }
})

test_that(".permute_treatment respects blocks (within-block shuffle)", {
  set.seed(1)
  blocks <- rep(1:5, each = 20)
  T_orig <- rep(0:1, each = 50)
  T_perm <- MLbalance:::.permute_treatment(T_orig, blocks = blocks)
  expect_equal(length(T_perm), length(T_orig))
  # Within each block, the set of treatment values is preserved
  for (b in unique(blocks)) {
    idx <- which(blocks == b)
    expect_equal(sort(T_perm[idx]), sort(T_orig[idx]))
  }
})

test_that(".permute_treatment handles clusters + blocks", {
  set.seed(1)
  n <- 100
  clusters <- rep(1:20, each = 5)
  blocks <- rep(1:2, each = 50)
  T_orig <- rep(rep(0:1, each = 10), each = 5)
  T_perm <- MLbalance:::.permute_treatment(T_orig, clusters = clusters, blocks = blocks)
  expect_equal(length(T_perm), n)
  # Constant within clusters
  for (cl in unique(clusters)) {
    vals <- T_perm[clusters == cl]
    expect_true(length(unique(vals)) == 1)
  }
})

# ============================================================================
# UTILS: .validate_clusters_blocks
# ============================================================================

test_that(".validate_clusters_blocks errors on wrong-length blocks", {
  expect_error(
    MLbalance:::.validate_clusters_blocks(NULL, rep(1, 5), n = 10, paired = FALSE),
    "same length"
  )
})

test_that(".validate_clusters_blocks errors on NA in blocks", {
  expect_error(
    MLbalance:::.validate_clusters_blocks(NULL, c(1, 2, NA, 4, 5), n = 5, paired = FALSE),
    "NA values"
  )
})

test_that(".validate_clusters_blocks errors on NA in clusters", {
  expect_error(
    MLbalance:::.validate_clusters_blocks(c(1, NA, 3), NULL, n = 3, paired = FALSE),
    "NA values"
  )
})

# ============================================================================
# FASTCPT: metrics and ensemble metrics
# ============================================================================

test_that("fastcpt works with metric = 'logscore'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", metric = "logscore",
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_equal(res$metric_name, "logscore")
})

test_that("fastcpt works with ensemble.metric = 'vote'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 ensemble.metric = "vote", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt works with ensemble.metric = 'mean.log'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 ensemble.metric = "mean.log", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt works with ensemble.metric = 'mean.mse'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 ensemble.metric = "mean.mse", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

# ============================================================================
# FASTCPT: leaveout = 0 (in-sample test statistic)
# ============================================================================

test_that("fastcpt works with leaveout = 0 (in-sample)", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", leaveout = 0,
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

# ============================================================================
# FASTCPT: paired permutation test
# ============================================================================

test_that("fastcpt works with paired = TRUE", {
  set.seed(1995)
  n <- 100
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", paired = TRUE,
                 perm.N = 30, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

# ============================================================================
# FASTCPT: data.frame Z with factor columns
# ============================================================================

test_that("fastcpt handles data.frame Z with factor columns", {
  set.seed(1995)
  n <- 200
  Z <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt handles data.frame Z with Date column", {
  set.seed(1995)
  n <- 200
  Z <- data.frame(
    x1 = rnorm(n),
    x2 = as.Date("2020-01-01") + sample(1:365, n, replace = TRUE)
  )
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

# ============================================================================
# FASTCPT: warnings for NA in Z
# ============================================================================

test_that("fastcpt warns on NA in Z", {
  set.seed(1995)
  n <- 200
  p <- 3
  Z <- matrix(rnorm(n * p), n, p)
  Z[1, 1] <- NA
  W <- rep(c(1, 2), each = n / 2)

  expect_warning(
    fastcpt(Z, W, class.methods = "forest", perm.N = 20, progress = FALSE),
    "NA or NaN"
  )
})

# ============================================================================
# FASTCPT: Inf in Z error
# ============================================================================

test_that("fastcpt errors on Inf in Z (matrix)", {
  Z <- matrix(rnorm(100), 50, 2)
  Z[1, 1] <- Inf
  W <- rep(c(1, 2), each = 25)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", perm.N = 10, progress = FALSE),
    "Inf"
  )
})

test_that("fastcpt errors on Inf in Z (data.frame)", {
  Z <- data.frame(a = c(Inf, rnorm(49)), b = rnorm(50))
  W <- rep(c(1, 2), each = 25)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", perm.N = 10, progress = FALSE),
    "Inf"
  )
})

# ============================================================================
# FASTCPT: unknown classifier errors
# ============================================================================

test_that("fastcpt errors on unknown classification method", {
  Z <- matrix(rnorm(100), 50, 2)
  W <- rep(c(1, 2), each = 25)

  expect_error(
    fastcpt(Z, W, class.methods = "xgboost", perm.N = 10, progress = FALSE),
    "Unknown classification method"
  )
})

# ============================================================================
# FASTCPT: classifier.args for forest
# ============================================================================

test_that("fastcpt respects classifier.args for forest (num.trees)", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 50,
                 progress = FALSE, classifier.args = list(num.trees = 100L))

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt respects classifier.args for ferns (depth)", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50,
                 progress = FALSE, classifier.args = list(ferns = 100L, depth = 3L))

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt respects classifier.args for glmnet2 (nfolds, alpha)", {
  skip_if_not_installed("glmnet")
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "glmnet2", perm.N = 50,
                 progress = FALSE, classifier.args = list(nfolds = 3L, alpha = 0.8))

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

# ============================================================================
# FASTCPT: multi-arm treatment (3+ classes)
# ============================================================================

test_that("fastcpt works with 3-class treatment", {
  set.seed(1995)
  n <- 150
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c("A", "B", "C"), each = 50)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

# ============================================================================
# FASTCPT.PLOT: multi-classifier plot/summary/print
# ============================================================================

test_that("plot.fastcpt works with multi-classifier result", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 perm.N = 50, progress = FALSE)
  pl <- plot(res)

  expect_s3_class(pl, "ggplot")
})

test_that("summary.fastcpt works with multi-classifier result", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 perm.N = 50, progress = FALSE)

  expect_output(summary(res), "Individual Classifiers")
  expect_output(summary(res), "Ensemble")
})

test_that("print.fastcpt includes design info with clusters and blocks", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  block_id <- rep(1:5, each = n / 5)
  W <- rep(rep(c(1, 2), each = n_clusters / 2), each = n / n_clusters)
  Z <- matrix(rnorm(n * p), n, p)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
                 clusters = cluster_id, blocks = block_id)

  expect_output(print(res), "clusters")
  expect_output(print(res), "blocks")
})

test_that("summary.fastcpt displays cluster and block info", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(rep(c(1, 2), each = n_clusters / 2), each = n / n_clusters)
  Z <- matrix(rnorm(n * p), n, p)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
                 clusters = cluster_id)

  expect_output(summary(res), "Clusters")
})

# ============================================================================
# FASTCPT: alpha stored correctly and pass/fail logic
# ============================================================================

test_that("fastcpt stores alpha and it affects pass/fail display", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
                 alpha = 0.10)

  expect_equal(res$alpha, 0.10)
  # summary should display alpha = 0.10
  expect_output(summary(res), "0.10")
})

# ============================================================================
# BALANCE: print and summary for multi-arm (no Y)
# ============================================================================

test_that("print.balance works for multi-arm without Y", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Ctrl", "A", "B"), n, replace = TRUE)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y = NULL, W = W, X = X, control = "Ctrl", perm.N = 50)
  })

  expect_output(print(res), "treatment arms")
  expect_output(print(res), "No outcome Y provided")
})

test_that("summary.balance works for multi-arm with Y", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Ctrl", "A", "B"), n, replace = TRUE)
  Y <- (W == "A") * 0.3 + (W == "B") * 0.6 + rnorm(n)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y = Y, W = W, X = X, control = "Ctrl", perm.N = 50)
  })

  expect_output(summary(res), "COVARIATE BALANCE ASSESSMENT")
  expect_output(summary(res), "TREATMENT EFFECT ESTIMATION")
  expect_output(summary(res), "ESTIMATOR DIVERGENCE TESTS")
})

# ============================================================================
# BALANCE: plot multi-arm
# ============================================================================

test_that("plot.balance works for multi-arm with Y", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Ctrl", "A", "B"), n, replace = TRUE)
  Y <- (W == "A") * 0.3 + (W == "B") * 0.6 + rnorm(n)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y = Y, W = W, X = X, control = "Ctrl", perm.N = 50)
  })

  pl_ps <- plot(res, which = "pscores")
  expect_true(inherits(pl_ps, "ggplot"))

  pl_null <- plot(res, which = "null_dist")
  expect_true(inherits(pl_null, "ggplot"))

  pl_eff <- plot(res, which = "effects")
  expect_true(inherits(pl_eff, "ggplot"))

  pl_all <- plot(res, which = "all")
  expect_true(inherits(pl_all, "ggplot") || inherits(pl_all, "patchwork"))
})

test_that("plot.balance multi-arm without Y skips effects", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Ctrl", "A", "B"), n, replace = TRUE)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y = NULL, W = W, X = X, control = "Ctrl", perm.N = 50)
  })

  # "all" should resolve to just pscores + null_dist
  expect_warning(
    plot(res, which = "effects"),
    "Outcome Y not provided"
  )
})

# ============================================================================
# BALANCE: single plots (combined = FALSE)
# ============================================================================

test_that("plot.balance returns list when combined = FALSE", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + rnorm(n)

  res <- balance(Y, W, X, perm.N = 50)
  pl <- plot(res, which = "all", combined = FALSE)

  expect_type(pl, "list")
  expect_true(length(pl) >= 2)
})

# ============================================================================
# BALANCE: blocks with Y (covers block dummy augmentation for grf)
# ============================================================================

test_that("balance works with blocks and Y (block dummies in grf)", {
  set.seed(1995)
  n <- 250
  p <- 5
  block_id <- rep(1:5, each = n / 5)
  W <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)
  Y <- W * 0.5 + X[, 1] * 0.3 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50, blocks = block_id)

  expect_s3_class(res, "balance")
  expect_true(is.numeric(res$dim$estimate))
  expect_true(is.numeric(res$aipw$estimate))
  expect_equal(res$blocks, block_id)
})

# ============================================================================
# BALANCE: summary without Y (early return before estimates section)
# ============================================================================

test_that("summary.balance works without Y (no estimates section)", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50)

  out <- capture.output(summary(res))
  expect_true(any(grepl("COVARIATE BALANCE ASSESSMENT", out)))
  # Should NOT have treatment effect section
  expect_false(any(grepl("TREATMENT EFFECT ESTIMATION", out)))
})

# ============================================================================
# BALANCE: print with clusters and blocks
# ============================================================================

test_that("print.balance shows design info with clusters", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(rep(c(0, 1), each = n_clusters / 2), each = n / n_clusters)
  X <- matrix(rnorm(n * p), n, p)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50,
                 clusters = cluster_id)

  expect_output(print(res), "clusters")
})

# ============================================================================
# BALANCE: summary with clusters (mentions cluster-robust)
# ============================================================================

test_that("summary.balance mentions cluster-robust when clusters provided", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(rep(c(0, 1), each = n_clusters / 2), each = n / n_clusters)
  X <- matrix(rnorm(n * p), n, p)
  Y <- W * 0.5 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50, clusters = cluster_id)

  expect_output(summary(res), "Cluster-robust")
})

# ============================================================================
# BALANCE: class.method = "forest"
# ============================================================================

test_that("balance works with class.method = 'forest'", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50, class.method = "forest")

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
})

# ============================================================================
# BALANCE: multi-arm errors with glmnet2 or lm
# ============================================================================

test_that("balance errors on multi-arm with glmnet2", {
  set.seed(1995)
  n <- 150
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("A", "B", "C"), n, replace = TRUE)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    expect_error(
      balance(Y = NULL, W = W, X = X, perm.N = 50, class.method = "glmnet2"),
      "binary treatment"
    )
  })
})

test_that("balance errors on multi-arm with lm", {
  set.seed(1995)
  n <- 150
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("A", "B", "C"), n, replace = TRUE)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    expect_error(
      balance(Y = NULL, W = W, X = X, perm.N = 50, class.method = "lm"),
      "binary treatment"
    )
  })
})

# ============================================================================
# BALANCE: multi-arm auto-detects control (verbose message)
# ============================================================================

test_that("balance auto-detects control for multi-arm and messages", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Ctrl", "A", "B"), n, replace = TRUE)

  expect_message(
    balance(Y = NULL, W = W, X = X, perm.N = 50),
    "Multi-arm treatment detected"
  )
})

# ============================================================================
# BALANCE: explicit control with verbose message
# ============================================================================

test_that("balance messages when explicit control is set", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  expect_message(
    balance(Y = NULL, W = W, X = X, perm.N = 50, control = "0"),
    "Using '0' as control group"
  )
})

# ============================================================================
# RANDOM_CHECK: facet = TRUE
# ============================================================================

test_that("random_check works with facet = TRUE", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- random_check(W_real = W, X = X, facet = TRUE)
  })

  expect_type(res, "list")
  expect_s3_class(res$plot, "ggplot")
})

# ============================================================================
# RANDOM_CHECK: data.frame X with factors
# ============================================================================

test_that("random_check works with data.frame X containing factors", {
  set.seed(1995)
  n <- 250
  X <- data.frame(
    x1 = rnorm(n),
    x2 = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    x3 = ordered(sample(c("low", "mid", "high"), n, replace = TRUE),
                 levels = c("low", "mid", "high"))
  )
  W <- rbinom(n, 1, 0.5)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- random_check(W_real = W, X = X)
  })

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
})

# ============================================================================
# RANDOM_CHECK: verbose messages
# ============================================================================

test_that("random_check prints permutation message when no W_sim", {
  set.seed(1995)
  n <- 250
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  withr::with_options(list(MLbalance.verbose = TRUE), {
    expect_message(
      random_check(W_real = W, X = X),
      "Permutated Treatment Assignment"
    )
  })
})

test_that("random_check prints simulated message when W_sim provided", {
  set.seed(1995)
  n <- 250
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  W_sim <- rbinom(n, 1, 0.5)

  withr::with_options(list(MLbalance.verbose = TRUE), {
    expect_message(
      random_check(W_real = W, W_sim = W_sim, X = X),
      "Simulated Assignment Vector Provided"
    )
  })
})

# ============================================================================
# RANDOM_CHECK: input validation
# ============================================================================

test_that("random_check errors on NA in W_sim", {
  set.seed(1995)
  n <- 50
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  W_sim <- W
  W_sim[1] <- NA

  expect_error(
    random_check(W_real = W, W_sim = W_sim, X = X),
    "NA values"
  )
})

test_that("random_check warns on NA in X", {
  set.seed(1995)
  n <- 250
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  X[1, 1] <- NA
  W <- rbinom(n, 1, 0.5)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    expect_warning(
      random_check(W_real = W, X = X),
      "NA or NaN"
    )
  })
})

test_that("random_check errors on Inf in X", {
  set.seed(1995)
  X <- data.frame(a = c(Inf, rnorm(49)), b = rnorm(50))
  W <- rbinom(50, 1, 0.5)

  expect_error(
    random_check(W_real = W, X = X),
    "Inf"
  )
})

test_that("random_check errors when W_real is not a vector", {
  set.seed(1995)
  X <- matrix(rnorm(100), 50, 2)
  W <- matrix(rbinom(50, 1, 0.5), 25, 2)

  expect_error(random_check(W_real = W, X = X), "vector")
})

test_that("random_check errors when W_sim is not a vector", {
  set.seed(1995)
  X <- matrix(rnorm(100), 50, 2)
  W <- rbinom(50, 1, 0.5)
  W_sim <- matrix(rbinom(50, 1, 0.5), 25, 2)

  expect_error(random_check(W_real = W, W_sim = W_sim, X = X), "vector")
})

# ============================================================================
# RANDOM_CHECK: character treatment (coercion to numeric)
# ============================================================================

test_that("random_check handles non-numeric W_real", {
  set.seed(1995)
  n <- 250
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- ifelse(rbinom(n, 1, 0.5) == 1, "treat", "ctrl")

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- random_check(W_real = W, X = X)
  })

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
})

# ============================================================================
# BALANCE: factor expansion warning for high-cardinality
# ============================================================================

test_that("balance warns on high factor expansion ratio", {
  set.seed(1995)
  n <- 250
  # Create a single factor with many levels to trigger the >5x expansion warning
  X <- data.frame(
    x1 = factor(sample(paste0("lev", 1:30), n, replace = TRUE))
  )
  W <- rbinom(n, 1, 0.5)

  # model.matrix creates 30 columns from 1, ratio = 30 > 5
  expect_warning(
    balance(Y = NULL, W = W, X = X, perm.N = 50),
    "Factor expansion"
  )
})

# ============================================================================
# BALANCE: Date/POSIXt columns in X
# ============================================================================

test_that("balance converts Date columns in X to numeric", {
  set.seed(1995)
  n <- 250
  X <- data.frame(
    x1 = rnorm(n),
    x2 = as.Date("2020-01-01") + sample(1:365, n, replace = TRUE)
  )
  W <- rbinom(n, 1, 0.5)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50)

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
})

# ============================================================================
# BALANCE: imp.predictors is returned
# ============================================================================

test_that("balance returns variable importance (imp.predictors)", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", 1:p)
  W <- rbinom(n, 1, 0.5)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50)

  expect_s3_class(res$imp.predictors, "data.frame")
  expect_true("varname" %in% names(res$imp.predictors))
  expect_true("vip" %in% names(res$imp.predictors))
})

# ============================================================================
# BALANCE: n, n_treated, n_control fields
# ============================================================================

test_that("balance returns correct sample size fields", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- c(rep(0, 100), rep(1, 150))

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50)

  expect_equal(res$n, 250)
  expect_equal(res$n_treated, 150)
  expect_equal(res$n_control, 100)
})

# ============================================================================
# BALANCE: multi-arm n_per_arm field
# ============================================================================

test_that("balance returns n_per_arm for multi-arm", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Ctrl", "A", "B"), n, replace = TRUE)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y = NULL, W = W, X = X, control = "Ctrl", perm.N = 50)
  })

  expect_true(!is.null(res$n_per_arm))
  expect_true("A" %in% names(res$n_per_arm))
  expect_true("B" %in% names(res$n_per_arm))
})

# ============================================================================
# FASTCPT: .softmax internal (via ferns in-sample)
# ============================================================================

test_that("fastcpt with leaveout = 0 exercises in-sample .softmax path", {
  set.seed(1995)
  n <- 100
  p <- 3
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  # leaveout = 0 uses testistrain = TRUE, which calls .softmax(fern.fit) (no newdata)
  res <- fastcpt(Z, W, class.methods = "ferns", leaveout = 0,
                 perm.N = 30, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$teststat))
})
