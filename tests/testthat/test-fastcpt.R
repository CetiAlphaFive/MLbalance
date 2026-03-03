test_that("fastcpt runs and returns expected structure", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_true(is.numeric(res$teststat))
  expect_true(is.numeric(res$nulldist))
  expect_length(res$nulldist, 50)
})

test_that("plot.fastcpt returns a ggplot", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE)
  pl <- plot(res)

  expect_s3_class(pl, "ggplot")
})

test_that("print.fastcpt runs without error", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE)

  expect_output(print(res), "Classification Permutation Test")
})

test_that("summary.fastcpt runs without error", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE)

  expect_output(summary(res), "Classification Permutation Test Results")
})
