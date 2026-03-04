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
