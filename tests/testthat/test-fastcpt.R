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

test_that("fastcpt works with clusters", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  W <- rep(c(1, 2), each = n / 2)
  # Ensure treatment is constant within clusters
  W <- rep(rep(c(1, 2), each = n_clusters / 2), each = n / n_clusters)
  Z <- matrix(rnorm(n * p), n, p)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
                 clusters = cluster_id)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_equal(res$clusters, cluster_id)
})

test_that("fastcpt works with blocks", {
  set.seed(1995)
  n <- 250
  p <- 5
  block_id <- rep(1:5, each = n / 5)
  W <- rep(c(1, 2), each = n / 2)
  Z <- matrix(rnorm(n * p), n, p)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
                 blocks = block_id)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_equal(res$blocks, block_id)
})

test_that("fastcpt works with both clusters and blocks", {
  set.seed(1995)
  n <- 250
  p <- 5
  n_clusters <- 50
  cluster_id <- rep(1:n_clusters, each = n / n_clusters)
  block_id <- rep(1:5, each = n / 5)
  # Treatment constant within clusters
  W <- rep(rep(c(1, 2), each = n_clusters / 2), each = n / n_clusters)
  Z <- matrix(rnorm(n * p), n, p)

  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
                 clusters = cluster_id, blocks = block_id)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt errors on paired + blocks", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)
  block_id <- rep(1:5, each = n / 5)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
            paired = TRUE, blocks = block_id),
    "Cannot use both"
  )
})

test_that("fastcpt errors on wrong-length clusters", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
            clusters = rep(1, 10)),
    "same length"
  )
})

test_that("fastcpt errors on treatment varying within cluster", {
  set.seed(1995)
  n <- 250
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)
  # Clusters that straddle treatment groups
  cluster_id <- rep(1:(n / 2), each = 2)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", perm.N = 50, progress = FALSE,
            clusters = cluster_id),
    "not constant within cluster"
  )
})

# --- New tests below ---

test_that("fastcpt works with forest classifier", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "forest", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("fastcpt works with glmnet2 classifier", {
  skip_if_not_installed("glmnet")
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "glmnet2", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("fastcpt works with multiple classifiers (ensemble)", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
  expect_true(is.matrix(res$nulldist))
  expect_equal(ncol(res$nulldist), 3)  # ferns, forest, ensemble
  expect_true("ensemble" %in% names(res$pvals))
})

test_that("fastcpt works with metric = 'rate'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", metric = "rate",
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt works with metric = 'mse'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "ferns", metric = "mse",
                 perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt works with comb.method = 'min'", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = c("ferns", "forest"),
                 comb.method = "min", perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt works with leaveout = 1 (LOO)", {
  set.seed(1995)
  n <- 50  # small for LOO speed
  p <- 3
  Z <- data.frame(matrix(rnorm(n * p), n, p))
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "forest", leaveout = 1,
                 leaveout.N = n, perm.N = 20, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt works with fractional leaveout", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- data.frame(matrix(rnorm(n * p), n, p))
  W <- rep(c(1, 2), each = n / 2)

  res <- fastcpt(Z, W, class.methods = "forest", leaveout = 0.2,
                 leaveout.N = 10, perm.N = 50, progress = FALSE)

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})

test_that("fastcpt errors on unknown metric", {
  set.seed(1995)
  Z <- matrix(rnorm(100), 50, 2)
  W <- rep(c(1, 2), each = 25)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", metric = "bogus",
            perm.N = 10, progress = FALSE),
    "Unknown metric"
  )
})

test_that("fastcpt errors on unknown ensemble metric", {
  set.seed(1995)
  Z <- matrix(rnorm(100), 50, 2)
  W <- rep(c(1, 2), each = 25)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", ensemble.metric = "bogus",
            perm.N = 10, progress = FALSE),
    "Unknown ensemble metric"
  )
})

test_that("fastcpt errors on unknown combination method", {
  set.seed(1995)
  Z <- matrix(rnorm(100), 50, 2)
  W <- rep(c(1, 2), each = 25)

  expect_error(
    fastcpt(Z, W, class.methods = "ferns", comb.method = "bogus",
            perm.N = 10, progress = FALSE),
    "Unknown combination method"
  )
})

test_that("fastcpt respects classifier.args for ferns", {
  set.seed(1995)
  n <- 200
  p <- 5
  Z <- matrix(rnorm(n * p), n, p)
  W <- rep(c(1, 2), each = n / 2)

  # Should run without error with custom fern count
  res <- fastcpt(Z, W, class.methods = "ferns", perm.N = 50,
                 progress = FALSE, classifier.args = list(ferns = 100L))

  expect_s3_class(res, "fastcpt")
  expect_true(is.numeric(res$pval))
})
