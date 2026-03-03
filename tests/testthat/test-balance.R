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
