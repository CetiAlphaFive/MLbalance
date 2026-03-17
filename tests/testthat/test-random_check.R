test_that("random_check runs and returns expected structure", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  res <- random_check(W_real = W, X = X)

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
  expect_length(res$treat.props.real, n)
  expect_true(is.numeric(res$treat.props.sim))
  expect_s3_class(res$plot, "ggplot")
})

test_that("vip returns a data frame with correct columns", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  Y <- rnorm(n)

  model <- grf::regression_forest(X = X, Y = Y, num.trees = 500, seed = 1995)
  vi <- vip(model)

  expect_s3_class(vi, "data.frame")
  expect_true("varname" %in% names(vi))
  expect_true("vip" %in% names(vi))
  expect_equal(nrow(vi), p)
})

test_that("random_check runs with clusters", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(rep(c(0, 1), each = n_clusters / 2), each = n / n_clusters)
  X <- matrix(rnorm(n * p), n, p)

  res <- random_check(W_real = W, X = X, clusters = cluster_id)

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
  expect_equal(res$clusters, cluster_id)
})

test_that("random_check runs with blocks", {
  set.seed(1995)
  n <- 250
  p <- 5
  block_id <- rep(1:5, each = n / 5)
  W <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)

  res <- random_check(W_real = W, X = X, blocks = block_id)

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
  expect_equal(res$blocks, block_id)
})

# --- New tests below ---

test_that("random_check works with explicit W_sim", {
  set.seed(1995)
  n <- 250
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W_real <- rbinom(n, 1, 0.5)
  W_sim <- rbinom(n, 1, 0.5)

  withr::with_options(list(MLbalance.verbose = FALSE), {
    res <- random_check(W_real = W_real, W_sim = W_sim, X = X)
  })

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
  expect_true(is.numeric(res$treat.props.sim))
  expect_s3_class(res$plot, "ggplot")
})

test_that("random_check works with clusters and blocks together", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  block_id <- rep(1:5, each = n / 5)
  W <- rep(rep(c(0, 1), each = n_clusters / 2), each = n / n_clusters)
  X <- matrix(rnorm(n * p), n, p)

  res <- random_check(W_real = W, X = X, clusters = cluster_id, blocks = block_id)

  expect_type(res, "list")
  expect_true(is.numeric(res$treat.props.real))
  expect_equal(res$clusters, cluster_id)
  expect_equal(res$blocks, block_id)
})
