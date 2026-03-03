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
