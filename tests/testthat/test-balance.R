test_that("balance runs without outcome Y", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50)

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
  expect_null(res$dim)
  expect_null(res$aipw)
  expect_equal(res$alpha, 0.05)
})

test_that("balance runs with outcome Y", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + X[, 1] * 0.3 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50)

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
  expect_true(is.numeric(res$dim$estimate))
  expect_true(is.numeric(res$aipw$estimate))
  expect_false(res$multiarm)
})

test_that("print.balance runs without error", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50)

  expect_output(print(res), "Balance Assessment")
})

test_that("summary.balance runs without error", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50)

  expect_output(summary(res), "COVARIATE BALANCE ASSESSMENT")
})

test_that("plot.balance returns a ggplot", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50)
  pl <- plot(res, which = "pscores")

  expect_s3_class(pl, "ggplot")
})

test_that("balance runs with clusters (no Y)", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(c(0, 1), each = n / 2)
  # Ensure treatment constant within clusters
  W <- rep(rep(c(0, 1), each = n_clusters / 2), each = n / n_clusters)
  X <- matrix(rnorm(n * p), n, p)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50, clusters = cluster_id)

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
  expect_equal(res$clusters, cluster_id)
})

test_that("balance runs with clusters (with Y)", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(rep(c(0, 1), each = n_clusters / 2), each = n / n_clusters)
  X <- matrix(rnorm(n * p), n, p)
  Y <- W * 0.5 + X[, 1] * 0.3 + rnorm(n)

  res <- balance(Y = Y, W = W, X = X, perm.N = 50, clusters = cluster_id)

  expect_s3_class(res, "balance")
  expect_true(is.numeric(res$dim$estimate))
  expect_true(is.numeric(res$aipw$estimate))
  expect_equal(res$clusters, cluster_id)
})

test_that("balance runs with blocks", {
  set.seed(1995)
  n <- 250
  p <- 5
  block_id <- rep(1:5, each = n / 5)
  W <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)

  res <- balance(Y = NULL, W = W, X = X, perm.N = 50, blocks = block_id)

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
  expect_equal(res$blocks, block_id)
})

# --- New tests below ---

test_that("balance works with multi-arm treatment", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Control", "A", "B"), n, replace = TRUE)
  Y <- (W == "A") * 0.3 + (W == "B") * 0.6 + rnorm(n)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y, W, X, control = "Control", perm.N = 50)
  })

  expect_s3_class(res, "balance")
  expect_true(res$multiarm)
  expect_equal(res$control, "Control")
  expect_true(is.logical(res$passed))
  expect_true(all(c("A", "B") %in% res$arms))

  # Verify named list structure for per-arm estimates

  expect_true(is.list(res$dim))
  expect_true("A" %in% names(res$dim))
  expect_true("B" %in% names(res$dim))
  expect_true(is.numeric(res$dim[["A"]]$estimate))
  expect_true(is.numeric(res$dim[["B"]]$estimate))
})

test_that("balance respects explicit control argument", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- sample(c("Placebo", "DrugA", "DrugB"), n, replace = TRUE)
  Y <- rnorm(n)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- balance(Y, W, X, control = "Placebo", perm.N = 50)
  })

  expect_equal(res$control, "Placebo")
  expect_false("Placebo" %in% res$arms)
})

test_that("balance handles factor covariates (ordered + unordered)", {
  set.seed(1995)
  n <- 250
  X <- data.frame(
    x1 = rnorm(n),
    x2 = factor(sample(c("low", "mid", "high"), n, replace = TRUE),
                levels = c("low", "mid", "high"), ordered = TRUE),
    x3 = factor(sample(c("cat", "dog"), n, replace = TRUE))
  )
  W <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)

  res <- balance(Y, W, X, perm.N = 50)

  expect_s3_class(res, "balance")
  expect_true(is.numeric(res$dim$estimate))
})

test_that("balance plot returns ggplot for different which values", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + rnorm(n)

  res <- balance(Y, W, X, perm.N = 50)

  pl_null <- plot(res, which = "null_dist")
  expect_s3_class(pl_null, "ggplot")

  pl_effects <- plot(res, which = "effects")
  expect_s3_class(pl_effects, "ggplot")

  pl_all <- plot(res, which = "all")
  # combined plot returns patchwork
  expect_true(inherits(pl_all, "ggplot") || inherits(pl_all, "patchwork"))
})

test_that("balance known-DGP correctness check", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  true_ate <- 5
  Y <- W * true_ate + X[, 1] * 0.3 + rnorm(n)

  res <- balance(Y, W, X, perm.N = 50)

  # All estimators should contain truth within estimate +/- 3*SE
  for (est_name in c("dim", "ipw", "aipw", "oadj")) {
    est <- res[[est_name]]
    lower <- est$estimate - 3 * est$std.err
    upper <- est$estimate + 3 * est$std.err
    expect_true(true_ate >= lower && true_ate <= upper,
                info = sprintf("%s: truth %.1f not in [%.2f, %.2f]",
                               est_name, true_ate, lower, upper))
  }
})

test_that("balance forwards fastcpt.args correctly", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  # Pass custom leaveout through fastcpt.args
  res <- balance(Y = NULL, W = W, X = X, perm.N = 50,
                 fastcpt.args = list(leaveout = 0.2, leaveout.N = 10))

  expect_s3_class(res, "balance")
  expect_true(is.logical(res$passed))
})

test_that("balance num.trees argument is respected", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + rnorm(n)

  # Should run with custom num.trees without error
  res <- balance(Y, W, X, perm.N = 50, num.trees = 500)

  expect_s3_class(res, "balance")
  expect_true(is.numeric(res$aipw$estimate))
})

test_that("balance ate_cov is positive semi-definite", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  Y <- W * 0.5 + X[, 1] * 0.3 + rnorm(n)

  res <- balance(Y, W, X, perm.N = 50)

  # Covariance matrix should be PSD
  eigenvalues <- eigen(res$ate_cov)$values
  expect_true(all(eigenvalues >= -1e-10),
              info = sprintf("Min eigenvalue: %e", min(eigenvalues)))
})
