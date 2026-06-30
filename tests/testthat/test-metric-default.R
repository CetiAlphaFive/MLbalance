# Auto-selected metric: OOB-capable backends (forest, ferns) at leaveout = 0
# default to "rate" (OOB classification accuracy), matching cpt's OOB-rate mode.
# Everything else (test-set splits, non-OOB backends, explicit metric) stays
# "probability".

make_Z <- function(n = 200, p = 5) {
  set.seed(1995)
  matrix(stats::rnorm(n * p), n, p)
}
binary_W <- function(n = 200) rep(c(0, 1), each = n / 2)

test_that("forest at leaveout = 0 defaults to rate", {
  res <- fastcpt(make_Z(), binary_W(), class.methods = "forest",
                 perm.N = 30, progress = FALSE)
  expect_identical(res$metric_name, "rate")
})

test_that("ferns at leaveout = 0 defaults to rate", {
  res <- fastcpt(make_Z(), binary_W(), class.methods = "ferns",
                 perm.N = 30, progress = FALSE)
  expect_identical(res$metric_name, "rate")
})

test_that("forest with leaveout > 0 stays probability", {
  res <- fastcpt(make_Z(), binary_W(), class.methods = "forest",
                 leaveout = 0.2, leaveout.N = 3, perm.N = 10, progress = FALSE)
  expect_identical(res$metric_name, "probability")
})

test_that("explicit metric is always respected", {
  res <- fastcpt(make_Z(), binary_W(), class.methods = "forest",
                 metric = "probability", perm.N = 30, progress = FALSE)
  expect_identical(res$metric_name, "probability")
})

test_that("non-OOB backend at leaveout = 0 stays probability", {
  res <- fastcpt(make_Z(), binary_W(), class.methods = "lm",
                 perm.N = 30, progress = FALSE)
  expect_identical(res$metric_name, "probability")
})

test_that("mixed OOB + non-OOB backends stay probability", {
  res <- fastcpt(make_Z(), binary_W(), class.methods = c("ferns", "lm"),
                 perm.N = 30, progress = FALSE)
  expect_identical(res$metric_name, "probability")
})

test_that("balance() with forest reports rate as the OOB statistic", {
  bal <- balance(W = binary_W(), X = make_Z(), class.method = "forest",
                 perm.N = 30)
  expect_identical(bal$balance_test$metric_name, "rate")
})
