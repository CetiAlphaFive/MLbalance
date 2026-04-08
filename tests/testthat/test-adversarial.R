# test-adversarial.R
# Adversarial tests: things a social scientist would actually do to break this package.
# Covers common real-world data mishaps, argument confusion, and edge cases.
#
# BUGS FOUND (tests that currently FAIL):
#
# BUG #1 — Character columns in covariates crash rFerns
#   balance() converts chars→factors for grf (via .prepare_covariates_for_grf)
#   but passes raw X_df with character columns to fastcpt(), where rFerns dies:
#   "All attributes must be either numeric or factor."
#   FIX: convert character→factor in X_df before calling fastcpt(), or in fastcpt() itself.
#   Affects: 4 tests (char covariates, mixed types, tibble+char, fastcpt+char)
#
# BUG #2 — High-cardinality factors (>30 levels) crash rFerns
#   rFerns hard limit: "unordered factor(s) with above 30 levels."
#   FIX: convert high-cardinality factors to one-hot before ferns, or error with
#   an informative message upstream (e.g., "Collapse factors with >30 levels").
#   Affects: 1 test (50-state factor)
#
# BUG #3 — Inf/-Inf in X not caught by validation
#   is.na(Inf) returns FALSE, so Inf slips through the NA/NaN check.
#   FIX: add `if (any(!is.finite(as.matrix(X[vapply(X, is.numeric, logical(1))]))))`
#   or simpler: check `any(is.infinite(as.matrix(X)))` alongside is.na().
#   Affects: 3 tests (Inf, -Inf, 1/x scenario)
#
# BUG #4 — Inf/-Inf in Y not caught by validation
#   Same as #3 but for outcome. Causes grf tuning crash:
#   "missing value where TRUE/FALSE needed" (sd(Inf) → NaN).
#   FIX: add is.infinite() check to Y validation.
#   Affects: 2 tests (Inf in Y, -Inf in Y)
#
# BUG #5 — Small n or extreme imbalance crashes grf honesty splitting
#   grf::boosted_regression_forest with default honesty = TRUE crashes when
#   effective sample is too small: "honesty fraction is too close to 1 or 0".
#   FIX: catch this and either (a) set honesty = FALSE for small samples,
#   (b) error with minimum-n guidance, or (c) wrap in tryCatch.
#   Affects: 3 tests (n=20, 190:10 imbalance, p>n with n=40)
#
# BUG #6 — random_check() crashes on character treatment
#   Passes character W_real directly to grf::boosted_regression_forest, which
#   requires numeric: "Observations must be numeric. GRF does not currently
#   support non-numeric observations."
#   FIX: convert W_real to numeric (0/1) inside random_check() before grf call.
#   Affects: 1 test (random_check with character treatment)
#
# BUG #7 — vip() crashes on boosted_regression_forest models
#   grf::variable_importance() → split_frequencies() hits:
#   "Index out of bounds: [index='_ci_group_size']"
#   Likely because variable_importance() doesn't support the boosted wrapper.
#   FIX: extract the underlying forest or use a different VIP method for boosted models.
#   Affects: 1 test (vip structure test)
#
# BUG #1 (expanded) — Date columns in X also crash rFerns
#   Date objects are neither numeric nor factor. Same root cause as char columns.
#   Affects: 1 test (Date column in covariates)

# ============================================================================
# Shared helpers
# ============================================================================

make_clean_data <- function(n = 200, p = 5) {
  set.seed(1995)
  W <- rep(0:1, each = n / 2)
  X <- data.frame(matrix(rnorm(n * p), n, p))
  Y <- rnorm(n) + W * 0.5
  list(W = W, X = X, Y = Y)
}

# ============================================================================
# balance(): type coercion — things people actually pass in
# ============================================================================

test_that("balance() handles character treatment ('treated'/'control')", {
  d <- make_clean_data()
  d$W <- ifelse(d$W == 1, "treated", "control")
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles logical treatment (TRUE/FALSE)", {
  d <- make_clean_data()
  d$W <- as.logical(d$W)
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles factor treatment", {
  d <- make_clean_data()
  d$W <- factor(d$W, labels = c("Control Group", "Treatment Group"))
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles numeric treatment coded as 1/2 instead of 0/1", {
  d <- make_clean_data()
  d$W <- d$W + 1L  # now 1 and 2

  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles tibble for X", {
  skip_if_not_installed("tibble")
  d <- make_clean_data()
  d$X <- tibble::as_tibble(d$X)
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles tibble for X with character columns", {
  # BUG #1: character columns in X crash rFerns — balance() converts chars→factors
  # for grf via .prepare_covariates_for_grf() but passes raw X_df to fastcpt(),
  # where rFerns chokes: "All attributes must be either numeric or factor"
  skip_if_not_installed("tibble")
  d <- make_clean_data(n = 200, p = 3)
  d$X$gender <- sample(c("Male", "Female"), 200, replace = TRUE)
  d$X <- tibble::as_tibble(d$X)
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

# ============================================================================
# balance(): covariate disasters
# ============================================================================

test_that("balance() handles matrix X without column names", {
  d <- make_clean_data()
  X_mat <- as.matrix(d$X)
  colnames(X_mat) <- NULL
  result <- balance(W = d$W, X = X_mat, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() errors when X is a single vector (no nrow)", {
  d <- make_clean_data()
  # User passes a vector instead of data.frame — nrow() is NULL
  expect_error(balance(W = d$W, X = d$X[, 1]))
})

test_that("balance() handles single-column data.frame for X", {
  d <- make_clean_data()
  result <- balance(W = d$W, X = d$X[, 1, drop = FALSE], perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles character covariates (state names, gender)", {
  # BUG #1: same char→ferns crash. balance() needs to convert chars→factors
  # in X_df before passing to fastcpt()
  d <- make_clean_data(n = 200, p = 2)
  d$X$state <- sample(c("CA", "TX", "NY", "FL", "OH"), 200, replace = TRUE)
  d$X$gender <- sample(c("Male", "Female", "Non-binary"), 200, replace = TRUE)
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles ordered factor covariate (education level)", {
  d <- make_clean_data(n = 200, p = 2)
  d$X$education <- ordered(
    sample(c("HS", "BA", "MA", "PhD"), 200, replace = TRUE),
    levels = c("HS", "BA", "MA", "PhD")
  )
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles mixed numeric, character, and factor covariates", {
  # BUG #1: same char→ferns crash. The "polsci survey data" scenario —
  # every grad student's data looks like this
  set.seed(1995)
  n <- 200
  W <- rep(0:1, each = 100)
  X <- data.frame(
    age           = rnorm(n, 45, 10),
    income        = rnorm(n, 50000, 15000),
    party         = sample(c("Dem", "Rep", "Ind"), n, replace = TRUE),
    education     = ordered(sample(c("HS", "BA", "MA", "PhD"), n, replace = TRUE),
                            levels = c("HS", "BA", "MA", "PhD")),
    female        = sample(0:1, n, replace = TRUE),
    state         = sample(state.abb[1:10], n, replace = TRUE)
  )
  result <- balance(W = W, X = X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() warns or handles constant covariate column", {
  d <- make_clean_data(n = 200, p = 3)
  d$X$constant_col <- 1  # column with zero variance
  # Should at least not crash
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles duplicate column names in X", {
  d <- make_clean_data(n = 200, p = 4)
  names(d$X) <- c("age", "age", "income", "income")
  # This is a real thing people do with copy-paste cbind()
  result <- balance(W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() errors on high-cardinality factor (50 states) with ferns", {
  # BUG #2: rFerns hard limit on unordered factors with >30 levels.
  # Now gives an informative error upstream.
  d <- make_clean_data(n = 200, p = 2)
  d$X$state <- factor(sample(state.abb, 200, replace = TRUE))
  expect_error(balance(W = d$W, X = d$X, perm.N = 50), "30 levels")
})

# ============================================================================
# balance(): NA / NaN / Inf — the silent data killers
# ============================================================================

test_that("balance() errors on NA in W", {
  d <- make_clean_data()
  d$W[5] <- NA
  expect_error(balance(W = d$W, X = d$X, perm.N = 50), "NA")
})

test_that("balance() warns on NA in X", {
  d <- make_clean_data()
  d$X[3, 2] <- NA
  expect_warning(balance(W = d$W, X = d$X, perm.N = 50), "NA")
})

test_that("balance() warns on NaN in X", {
  d <- make_clean_data()
  d$X[1, 1] <- NaN
  expect_warning(balance(W = d$W, X = d$X, perm.N = 50), "NA|NaN")
})

test_that("balance() errors on Inf in X (from log(0) or 1/0)", {
  # BUG #3: is.na() doesn't catch Inf. Validation uses:
  #   if (any(is.na(as.matrix(X)))) stop(...)
  # but is.na(Inf) is FALSE. Inf slips through and causes cryptic
  # downstream failures. Need: is.infinite() or is.finite() check.
  d <- make_clean_data()
  d$X[1, 1] <- Inf
  expect_error(balance(W = d$W, X = d$X, perm.N = 50), "Inf|infinite")
})

test_that("balance() errors on -Inf in X", {
  # BUG #3: same Inf validation gap
  d <- make_clean_data()
  d$X[1, 1] <- -Inf
  expect_error(balance(W = d$W, X = d$X, perm.N = 50), "Inf|infinite")
})

test_that("balance() errors on NA in Y", {
  d <- make_clean_data()
  d$Y[10] <- NA
  expect_error(balance(Y = d$Y, W = d$W, X = d$X, perm.N = 50), "NA")
})

test_that("balance() errors on Inf in Y (from log(0))", {
  # BUG #4: Y validation checks is.na(Y) but not is.infinite(Y).
  # Inf in Y causes grf tuning to crash with:
  #   "missing value where TRUE/FALSE needed" (sd() returns NA on Inf)
  d <- make_clean_data()
  d$Y[1] <- Inf
  expect_error(balance(Y = d$Y, W = d$W, X = d$X, perm.N = 50), "Inf|infinite")
})

# ============================================================================
# balance(): wrong shapes and argument confusion
# ============================================================================

test_that("balance() errors when W and X have different lengths (off-by-one)", {
  d <- make_clean_data(n = 200)
  # User filtered X but forgot to filter W
  expect_error(
    balance(W = d$W, X = d$X[1:199, ]),
    "same number"
  )
})

test_that("balance() errors when Y and W have different lengths", {
  d <- make_clean_data(n = 200)
  expect_error(
    balance(Y = d$Y[1:150], W = d$W, X = d$X, perm.N = 50),
    "same number"
  )
})

test_that("balance() errors when Y is character", {
  d <- make_clean_data()
  d$Y <- ifelse(d$Y > 0, "positive", "negative")
  expect_error(balance(Y = d$Y, W = d$W, X = d$X, perm.N = 50), "numeric")
})

test_that("balance() errors when Y is a factor", {
  d <- make_clean_data()
  d$Y <- factor(ifelse(d$Y > 0, "high", "low"))
  expect_error(balance(Y = d$Y, W = d$W, X = d$X, perm.N = 50), "numeric")
})

test_that("balance() errors when W has only one unique value (forgot control group)", {
  d <- make_clean_data()
  d$W <- rep(1, 200)
  expect_error(balance(W = d$W, X = d$X, perm.N = 50), "at least 2")
})

test_that("balance() errors when control level doesn't exist in W", {
  d <- make_clean_data()
  expect_error(
    balance(W = d$W, X = d$X, perm.N = 50, control = "placebo"),
    "not found"
  )
})

test_that("balance() errors when X is a plain list (not data.frame)", {
  d <- make_clean_data()
  X_list <- list(a = rnorm(200), b = rnorm(200))
  expect_error(balance(W = d$W, X = X_list, perm.N = 50))
})

# ============================================================================
# balance(): small sample and degenerate cases
# ============================================================================

test_that("balance() errors on very small n (n = 20, 10 per group)", {
  # BUG #5: grf crashes with honesty fraction error for small samples.
  # Now gives an informative error.
  set.seed(1995)
  n <- 20
  W <- rep(0:1, each = 10)
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  expect_error(balance(W = W, X = X, perm.N = 50), "too small|too imbalanced")
})

test_that("balance() handles extreme treatment imbalance (190 vs 10)", {
  # grf handles this sample size fine; overlap diagnostics should flag it
  set.seed(1995)
  n <- 200
  W <- c(rep(0, 190), rep(1, 10))
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  result <- balance(W = W, X = X, perm.N = 50)
  expect_s3_class(result, "balance")
  expect_true(result$overlap_flag)
})

test_that("balance() errors on p > n (more covariates than observations)", {
  # BUG #5: same grf honesty crash for small n (40 obs, 20 per group)
  # Now gives an informative error.
  set.seed(1995)
  n <- 40
  p <- 60
  W <- rep(0:1, each = 20)
  X <- data.frame(matrix(rnorm(n * p), n, p))
  expect_error(balance(W = W, X = X, perm.N = 50), "too small|too imbalanced")
})

test_that("balance() handles all identical covariate rows (no variation)", {
  set.seed(1995)
  n <- 100
  W <- rep(0:1, each = 50)
  X <- data.frame(x1 = rep(1, n), x2 = rep(2, n), x3 = rep(3, n))
  # Classifier has nothing to learn; should still return results
  result <- balance(W = W, X = X, perm.N = 50)
  expect_s3_class(result, "balance")
  # Balance should pass since covariates are identical across groups
  expect_true(result$passed)
})

test_that("balance() handles binary outcome (Y = 0/1)", {
  d <- make_clean_data(n = 200)
  d$Y <- rbinom(200, 1, 0.5 + 0.1 * d$W)
  result <- balance(Y = d$Y, W = d$W, X = d$X, perm.N = 50)
  expect_s3_class(result, "balance")
})

# ============================================================================
# balance(): multi-arm headaches
# ============================================================================

test_that("balance() handles multi-arm with very unequal group sizes", {
  set.seed(1995)
  n <- 200
  W <- c(rep("control", 150), rep("low_dose", 30), rep("high_dose", 20))
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  result <- balance(W = W, X = X, perm.N = 50, control = "control")
  expect_s3_class(result, "balance")
  expect_true(result$multiarm)
})

test_that("balance() handles multi-arm when user forgets control argument", {
  set.seed(1995)
  n <- 150
  W <- rep(c("placebo", "drug_a", "drug_b"), each = 50)
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  # Default: first factor level as control (alphabetical → drug_a)
  # User probably wanted "placebo" — tests that it at least runs
  expect_message(
    result <- balance(W = W, X = X, perm.N = 50),
    "Multi-arm|control"
  )
  expect_s3_class(result, "balance")
})

test_that("balance() handles multi-arm with numeric treatment (0, 1, 2)", {
  set.seed(1995)
  n <- 150
  W <- rep(0:2, each = 50)
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  result <- balance(W = W, X = X, perm.N = 50, control = "0")
  expect_s3_class(result, "balance")
})

# ============================================================================
# balance(): clusters and blocks edge cases
# ============================================================================

test_that("balance() errors when clusters have NA", {
  d <- make_clean_data(n = 200)
  clusters <- rep(1:40, each = 5)
  clusters[3] <- NA
  expect_error(
    balance(W = d$W, X = d$X, clusters = clusters, perm.N = 50),
    "NA"
  )
})

test_that("balance() handles singleton clusters (each obs is own cluster)", {
  d <- make_clean_data(n = 100)
  clusters <- seq_len(100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], clusters = clusters, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() handles blocks where one block has only 1 treated", {
  set.seed(1995)
  n <- 100
  blocks <- rep(1:5, each = 20)
  # Block 1: only 1 treated, 19 control; others balanced
  W <- c(c(0, rep(0, 18), 1), rep(c(rep(0, 10), rep(1, 10)), 4))
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  result <- balance(W = W, X = X, blocks = blocks, perm.N = 50)
  expect_s3_class(result, "balance")
})

test_that("balance() errors when clusters and treatment are mismatched", {
  d <- make_clean_data(n = 100)
  # Clusters where treatment varies within cluster — invalid
  clusters <- rep(1:20, each = 5)
  W_bad <- sample(0:1, 100, replace = TRUE)
  expect_error(
    balance(W = W_bad, X = d$X[1:100, ], clusters = clusters, perm.N = 50),
    "not constant within"
  )
})

# ============================================================================
# balance(): argument foot-guns
# ============================================================================

test_that("balance() handles perm.N = 1 without crashing", {
  d <- make_clean_data(n = 100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], perm.N = 1)
  expect_s3_class(result, "balance")
})

test_that("balance() handles alpha = 0 (reject everything)", {
  d <- make_clean_data(n = 100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], perm.N = 50, alpha = 0)
  expect_s3_class(result, "balance")
  # alpha = 0 means p > 0 is always true, so should pass
  expect_true(result$passed)
})

test_that("balance() handles alpha = 1 (accept everything)", {
  d <- make_clean_data(n = 100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], perm.N = 50, alpha = 1)
  expect_s3_class(result, "balance")
})

test_that("balance() handles overlap.threshold = c(0, 1) (no flagging)", {
  d <- make_clean_data(n = 100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], perm.N = 50,
                    overlap.threshold = c(0, 1))
  expect_s3_class(result, "balance")
  expect_false(result$overlap_flag)
})

# ============================================================================
# balance(): S3 methods on weird results
# ============================================================================

test_that("print/summary/plot work on minimal balance object (no Y)", {
  d <- make_clean_data(n = 100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], perm.N = 50)
  expect_output(print(result))
  expect_output(summary(result))
  p <- plot(result, which = "null_dist")
  expect_s3_class(p, "gg")
})

test_that("plot.balance() handles which = 'effects' when Y is NULL", {
  d <- make_clean_data(n = 100)
  result <- balance(W = d$W[1:100], X = d$X[1:100, ], perm.N = 50)
  # No outcome → no effects to plot; currently warns (not errors)
  expect_warning(plot(result, which = "effects"), "Outcome Y not provided")
})

# ============================================================================
# fastcpt(): dumb inputs
# ============================================================================

test_that("fastcpt() errors when Z is a vector", {
  set.seed(1995)
  Z_vec <- rnorm(100)
  T_vec <- rep(0:1, each = 50)
  expect_error(fastcpt(Z = Z_vec, T = T_vec, perm.N = 50), "matrix or data frame")
})

test_that("fastcpt() errors on NA in treatment vector", {
  set.seed(1995)
  Z <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  T_vec <- rep(0:1, each = 50)
  T_vec[5] <- NA
  expect_error(fastcpt(Z = Z, T = T_vec, perm.N = 50), "NA")
})

test_that("fastcpt() handles T as character ('treat'/'ctrl')", {
  set.seed(1995)
  Z <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  T_vec <- rep(c("ctrl", "treat"), each = 50)
  result <- fastcpt(Z = Z, T = T_vec, perm.N = 50)
  expect_s3_class(result, "fastcpt")
})

test_that("fastcpt() handles logical T (TRUE/FALSE)", {
  set.seed(1995)
  Z <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  T_vec <- rep(c(TRUE, FALSE), each = 50)
  result <- fastcpt(Z = Z, T = T_vec, perm.N = 50)
  expect_s3_class(result, "fastcpt")
})

test_that("fastcpt() errors on bogus class.methods", {
  set.seed(1995)
  Z <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  T_vec <- rep(0:1, each = 50)
  expect_error(fastcpt(Z = Z, T = T_vec, class.methods = "xgboost", perm.N = 50))
})

test_that("fastcpt() handles Z with character columns (passed to ferns)", {
  # BUG #1: same character column issue — fastcpt() passes raw Z to ferns
  # which requires numeric or factor. Need char→factor conversion in fastcpt().
  set.seed(1995)
  n <- 100
  Z <- data.frame(
    x1 = rnorm(n),
    party = sample(c("Dem", "Rep", "Ind"), n, replace = TRUE)
  )
  T_vec <- rep(0:1, each = 50)
  result <- fastcpt(Z = Z, T = T_vec, perm.N = 50)
  expect_s3_class(result, "fastcpt")
})

test_that("fastcpt() handles single-column Z", {
  set.seed(1995)
  Z <- data.frame(x1 = rnorm(100))
  T_vec <- rep(0:1, each = 50)
  result <- fastcpt(Z = Z, T = T_vec, perm.N = 50)
  expect_s3_class(result, "fastcpt")
})

# ============================================================================
# random_check(): user confusion
# ============================================================================

test_that("random_check() errors when W_real has NAs", {
  set.seed(1995)
  W_real <- rep(0:1, each = 50)
  W_real[3] <- NA
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  expect_error(random_check(W_real = W_real, X = X), "NA")
})

test_that("random_check() errors when X is a vector", {
  set.seed(1995)
  W_real <- rep(0:1, each = 50)
  X_vec <- rnorm(100)
  expect_error(random_check(W_real = W_real, X = X_vec), "matrix or data frame")
})

test_that("random_check() handles character treatment", {
  # BUG #6: random_check() passes character W_real directly to grf which requires numeric.
  # Crashes: "Observations must be numeric."
  set.seed(1995)
  W_real <- rep(c("treated", "control"), each = 50)
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  result <- random_check(W_real = W_real, X = X)
  expect_type(result, "list")
  expect_s3_class(result$plot, "gg")
})

test_that("random_check() handles W_sim with different treatment coding", {
  set.seed(1995)
  n <- 100
  W_real <- rep(0:1, each = 50)
  W_sim <- sample(0:1, n, replace = TRUE)
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  result <- random_check(W_real = W_real, W_sim = W_sim, X = X)
  expect_type(result, "list")
})

# ============================================================================
# ============================================================================
# Cross-function: the "I just loaded my Stata data" scenario
# ============================================================================

test_that("balance() handles Stata-like data (haven-labelled would be numeric)", {
  # Simulating what read_dta() produces: numeric columns, character labels
  set.seed(1995)
  n <- 200
  W <- rep(0:1, each = 100)
  X <- data.frame(
    age = round(rnorm(n, 45, 12)),
    female = sample(0:1, n, replace = TRUE),
    educ = sample(1:5, n, replace = TRUE),      # coded as integers (1=HS, 2=BA, etc.)
    income_log = rnorm(n, 10.5, 0.8),
    pid7 = sample(1:7, n, replace = TRUE)        # 7-point party ID scale
  )
  Y <- rnorm(n)
  result <- balance(Y = Y, W = W, X = X, perm.N = 50)
  expect_s3_class(result, "balance")
  expect_false(is.null(result$dim))
  expect_false(is.null(result$aipw))
})

test_that("balance() handles data where user did log(0) creating -Inf in Y", {
  # BUG #4: same Inf validation gap in Y
  d <- make_clean_data()
  d$Y[5] <- -Inf
  expect_error(balance(Y = d$Y, W = d$W, X = d$X, perm.N = 50), "Inf|infinite")
})

test_that("balance() handles data where user did 1/x creating Inf in X", {
  # BUG #3: same Inf validation gap in X
  d <- make_clean_data()
  d$X[10, 1] <- Inf  # from 1/0
  expect_error(balance(W = d$W, X = d$X, perm.N = 50), "Inf|infinite")
})

# ============================================================================
# Date and POSIXct columns — left in by accident
# ============================================================================

test_that("balance() handles Date column in covariates", {
  # BUG #1 (variant): Date objects are neither numeric nor factor — rFerns crashes.
  # Common scenario: user forgets to drop survey_date before running balance().
  set.seed(1995)
  n <- 200
  W <- rep(0:1, each = 100)
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    survey_date = as.Date("2024-01-01") + sample(0:365, n, replace = TRUE)
  )
  result <- balance(W = W, X = X, perm.N = 50)
  expect_s3_class(result, "balance")
})
