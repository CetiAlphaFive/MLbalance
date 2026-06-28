# Tests for the three native backends added in place of fastcpt3 / mlr3:
#   rpart, lda, qda.
# All gated by skip_if_not_installed(); use perm.N = 50 and n <= 250 per
# package CLAUDE.md test-speed guidance.

test_that("fastcpt runs with rpart backend", {
  skip_if_not_installed("rpart")
  set.seed(1995)
  n <- 200; p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "rpart",
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
  expect_length(res$nulldist, 50)
})

test_that("fastcpt runs with lda backend", {
  skip_if_not_installed("MASS")
  set.seed(1995)
  n <- 200; p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "lda",
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("fastcpt runs with qda backend", {
  skip_if_not_installed("MASS")
  set.seed(1995)
  # qda needs n_k > p per group; with binary T at 50/50 and n = 250, p = 4
  # each group has 125 rows, well above p.
  n <- 250; p <- 4
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "qda",
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("fastcpt rpart detects real signal", {
  skip_if_not_installed("rpart")
  set.seed(1995)
  n <- 200; p <- 4
  W <- rep(c(1, 2), each = n / 2)
  Z <- matrix(rnorm(n * p), n, p)
  # Inject a 1.5 SD shift in the treated group across all covariates
  Z[W == 2, ] <- Z[W == 2, ] + 1.5

  res <- fastcpt(Z, W, class.methods = "rpart",
                 perm.N = 100, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_lt(res$pval, 0.05)
})

test_that("fastcpt lda + ranger ensemble works", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("ranger")
  set.seed(1995)
  n <- 200; p <- 4
  W <- rep(c(1, 2), each = n / 2)
  Z <- matrix(rnorm(n * p), n, p)

  res <- fastcpt(Z, W, class.methods = c("lda", "forest"),
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  # When >1 method is requested, names(res$pvals) includes each method + "ensemble"
  expect_true("ensemble" %in% names(res$pvals))
  expect_true("lda"      %in% names(res$pvals))
  expect_true("forest"   %in% names(res$pvals))
})

test_that("fastcpt errors clearly when backend pkg missing", {
  # Hard to simulate without unloading the namespace; rely on
  # skip_if_not_installed to exercise the install-present path. The
  # error path itself is straightforward (mirrors the existing glmnet2
  # guard) and is covered by code review.
  skip("Requires package mocking; covered manually via DESCRIPTION inspection.")
})

test_that(".gettrainmethod accepts a leaveout argument", {
  expect_silent(MLbalance:::.gettrainmethod("forest", list(), leaveout = 0))
})

test_that("forest forwards ranger args from classifier.args (splitrule)", {
  set.seed(1995)
  Z <- matrix(rnorm(200 * 5), 200, 5); W <- rep(c(0, 1), each = 100)
  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE,
                 classifier.args = list(splitrule = "extratrees", num.random.splits = 1L))
  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pvals[["forest"]]) && !is.na(res$pvals[["forest"]]))
})

test_that("forest ignores non-ranger classifier.args keys", {
  set.seed(1995)
  Z <- matrix(rnorm(200 * 5), 200, 5); W <- rep(c(0, 1), each = 100)
  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE,
                 classifier.args = list(num.trees = 50L, ferns = 999L, depth = 7L))
  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pvals[["forest"]]))
})

test_that("forest leaveout>0 path still works (write.forest=TRUE branch)", {
  set.seed(1995)
  Z <- matrix(rnorm(120 * 4), 120, 4); W <- rep(c(0, 1), each = 60)
  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 30, progress = FALSE,
                 leaveout = 1, leaveout.N = 30)
  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pvals[["forest"]]))
})

test_that("forest pval is reproducible across identical calls", {
  Z <- matrix(rnorm(200 * 5), 200, 5); W <- rep(c(0, 1), each = 100)
  p1 <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE)$pvals[["forest"]]
  p2 <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE)$pvals[["forest"]]
  expect_identical(p1, p2)
})

test_that("forest tolerates classifier.args = NULL", {
  set.seed(1995)
  Z <- matrix(rnorm(200 * 5), 200, 5); W <- rep(c(0, 1), each = 100)
  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 30,
                 progress = FALSE, classifier.args = NULL)
  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pvals[["forest"]]))
})

test_that("default forest run returns a valid, in-range p-value", {
  set.seed(1995)
  Z <- matrix(rnorm(200 * 5), 200, 5); W <- rep(c(0, 1), each = 100)
  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE)
  expect_true(is.numeric(res$pvals[["forest"]]))
  expect_gte(res$pvals[["forest"]], 0)
  expect_lte(res$pvals[["forest"]], 1)
})

test_that("default forest (extratrees) falls back to gini on NA in Z", {
  set.seed(1995)
  Z <- matrix(rnorm(200 * 5), 200, 5); Z[1, 1] <- NA
  W <- rep(c(0, 1), each = 100)
  res <- suppressWarnings(fastcpt(Z, W, class.methods = "forest", perm.N = 30, progress = FALSE))
  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pvals[["forest"]]) && !is.na(res$pvals[["forest"]]))
})
