# Getting Started with MLbalance

## Introduction

MLbalance provides machine learning-based covariate balance tests and
causal effect estimation for experimental and observational data. The
core tool is a **classification permutation test (CPT)** based on
Gagnon-Bartsch & Shem-Tov (2019): if a classifier can distinguish
treated from control units better than chance, covariate balance fails.

The package also estimates average treatment effects (ATEs) using four
approaches — difference-in-means, IPW, outcome-adjusted, and AIPW — all
implemented via
[`grf::causal_forest`](https://rdrr.io/pkg/grf/man/causal_forest.html).

## Basic Usage: Binary Treatment

The main function is
[`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md),
which runs the balance test and (optionally) estimates treatment effects
in a single call.

``` r
library(MLbalance)

set.seed(1995)
n <- 500
p <- 10
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
Y <- W * 0.5 + X[, 1] * 0.3 + rnorm(n)

result <- balance(Y, W, X, perm.N = 200)
result
summary(result)
```

The [`summary()`](https://rdrr.io/r/base/summary.html) output includes:

- The CPT result (pass/fail) with permutation null distribution
  statistics
- Propensity score diagnostics (real vs. null distributions)
- Four ATE estimates with standard errors and confidence intervals
- Estimator divergence tests that flag when different estimators
  disagree

### Plotting

``` r
plot(result)
```

The combined plot shows three panels: propensity score distributions,
the CPT null distribution, and treatment effect estimates with
confidence intervals.

You can also request individual panels:

``` r
plot(result, which = "pscores")
plot(result, which = "null_dist")
plot(result, which = "effects")
```

## Multi-Arm Designs

For experiments with more than two treatment arms,
[`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md)
runs a joint K-class CPT on the full data, then performs pairwise binary
comparisons against the control group for ATE estimation.

``` r
W_multi <- sample(c("Control", "Drug A", "Drug B"), n, replace = TRUE)
Y_multi <- (W_multi == "Drug A") * 0.3 +
           (W_multi == "Drug B") * 0.6 + rnorm(n)

result_multi <- balance(Y_multi, W_multi, X, control = "Control", perm.N = 200)
summary(result_multi)
plot(result_multi)
```

Use the `control` argument to specify which level serves as the
reference group. If omitted, the first factor level is used.

## Cluster Randomization and Blocking

When treatment is assigned at the cluster level, pass a `clusters`
vector. Permutations will shuffle treatment labels at the cluster level,
and standard errors will account for clustering.

``` r
n_clusters <- 100
cluster_id <- rep(seq_len(n_clusters), each = 5)
W_clust <- rep(rbinom(n_clusters, 1, 0.5), each = 5)
X_clust <- matrix(rnorm(500 * p), 500, p)
Y_clust <- W_clust * 0.5 + rnorm(500)

result_clust <- balance(Y_clust, W_clust, X_clust,
                        clusters = cluster_id, perm.N = 200)
```

For block-randomized designs, pass a `blocks` vector. Permutations are
restricted to within each block.

``` r
block_id <- rep(1:5, each = 100)
W_block <- rbinom(500, 1, 0.5)
X_block <- matrix(rnorm(500 * p), 500, p)
Y_block <- W_block * 0.5 + rnorm(500)

result_block <- balance(Y_block, W_block, X_block,
                        blocks = block_id, perm.N = 200)
```

Both `clusters` and `blocks` can be used simultaneously for
cluster-randomized block designs.

## Classifier Selection

The CPT supports three classifier backends, selected via the
`class.method` argument in
[`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md)
(or `class.methods` in
[`fastcpt()`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md)):

| Backend             | Speed     | Interactions     | Multi-class | Best for                    |
|---------------------|-----------|------------------|-------------|-----------------------------|
| `"ferns"` (default) | Very fast | Native           | Yes         | General use, large datasets |
| `"forest"`          | Moderate  | Learned          | Yes         | Strong nonlinearities       |
| `"glmnet2"`         | Fast      | 2-way (explicit) | Binary only | Linear/low-dim settings     |

The `"ferns"` backend (random ferns) is recommended for most
applications: it handles interactions natively, runs very quickly, and
works well with both binary and multi-arm treatments.

For advanced control, pass classifier hyperparameters via
`fastcpt.args`:

``` r
# Custom ranger forest with 1000 trees
result <- balance(Y, W, X, class.method = "forest",
                  fastcpt.args = list(classifier.args = list(num.trees = 1000)),
                  perm.N = 200)

# Custom glmnet with 10-fold CV
result <- balance(Y, W, X, class.method = "glmnet2",
                  fastcpt.args = list(classifier.args = list(nfolds = 10)),
                  perm.N = 200)
```

## Interpreting Overlap Diagnostics

When propensity scores are extreme (near 0 or 1), standard AIPW
estimators can become unstable.
[`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md)
automatically flags overlap issues and computes overlap-weighted (OW)
estimates as a fallback (Li, Morgan & Zaslavsky, 2018).

The overlap threshold is configurable:

``` r
result <- balance(Y, W, X, overlap.threshold = c(0.10, 0.90), perm.N = 200)
```

## Estimator Divergence Tests

The [`summary()`](https://rdrr.io/r/base/summary.html) output includes
pairwise divergence tests between the four ATE estimators. These test
whether the difference between two estimators is statistically
significant, using influence-function-based covariance estimation.

Divergence signals that the choice of nuisance model specification
matters for this dataset:

- **DiM vs IPW:** Propensity reweighting changes the estimate
- **DiM vs Outcome-adjusted:** Outcome regression adjustment changes the
  estimate
- **DiM vs AIPW:** Combined adjustment changes the estimate
- **IPW vs AIPW:** Outcome modeling adds information beyond propensity
  reweighting alone

Agreement across all four estimators indicates robustness to modeling
assumptions.

## Using fastcpt Directly

For a standalone balance test without treatment effect estimation, use
[`fastcpt()`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md)
directly:

``` r
Z <- matrix(rnorm(500 * 10), 500, 10)
T_var <- rep(c(1, 2), each = 250)

cpt_result <- fastcpt(Z, T_var, class.methods = "ferns", perm.N = 500)
summary(cpt_result)
plot(cpt_result)
```

[`fastcpt()`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md)
supports multiple classifiers simultaneously with ensemble combination:

``` r
cpt_ensemble <- fastcpt(Z, T_var,
                        class.methods = c("ferns", "forest"),
                        perm.N = 500)
summary(cpt_ensemble)
```
