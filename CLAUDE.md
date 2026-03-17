# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What This Package Does

MLbalance (v0.2) provides ML-based covariate balance tests and causal
effect estimation for experimental and observational data. The core tool
is a fast classification permutation test (CPT) based on Gagnon-Bartsch
& Shem-Tov (2019). If a classifier can distinguish treated from control
units better than chance, balance fails. The package also estimates ATEs
using four approaches: difference-in-means (DiM), inverse propensity
weighted (IPW), outcome-adjusted, and AIPW (doubly robust) (all via
`grf`).

## Build and Test Commands

``` bash
# Check package (full R CMD check)
R CMD build . && R CMD check MLbalance_*.tar.gz

# Or via devtools in R:
devtools::check()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-fastcpt.R")

# Document (regenerate NAMESPACE, Rd files)
devtools::document()

# Load package for interactive development
devtools::load_all()
```

Tests are slow (each
[`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md)
test fits multiple boosted forests). Use `perm.N = 50` and `n <= 250` in
tests.

## Architecture

Three-layer design, each layer an S3 class:

1.  **`fastcpt(Z, T, ...)`** (`R/fastcpt.R`) — The permutation test
    engine. Trains classifiers (ferns, ranger forest, or glmnet elastic
    net) on real treatment labels, then on `perm.N` permuted labels to
    build a null distribution. Returns p-value, test statistic, and null
    distribution. Supports `parallel = TRUE` via `mirai`. S3 methods in
    `R/fastcpt.plot.R`.

2.  **`balance(Y, W, X, ...)`** (`R/balance.R`) — The main user-facing
    function. Calls
    [`fastcpt()`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md)
    for the balance test, then fits
    [`grf::boosted_regression_forest`](https://rdrr.io/pkg/grf/man/boosted_regression_forest.html)
    for propensity/outcome models and
    [`grf::causal_forest`](https://rdrr.io/pkg/grf/man/causal_forest.html)
    for ATE estimation (DiM, IPW, outcome-adjusted, AIPW). Handles
    multi-arm treatments via pairwise comparisons against a control
    level, with a joint K-class CPT. Detects extreme propensity scores
    and provides overlap-weighted (OW) fallback estimates. Returns S3
    class “balance” with print/summary/plot methods.

3.  **`random_check(W_real, X, ...)`** (`R/random_check.R`) —
    Lightweight diagnostic. Fits boosted RF propensity models on real
    vs. permuted/simulated treatment, returns overlapping propensity
    score distributions. Also exports
    [`vip()`](https://cetialphafive.github.io/MLbalance/reference/vip.md)
    for variable importance from grf models.

**Shared utilities** (`R/utils.R`): - `.save_rng_state()` /
`.restore_rng_state()` — All exported functions save and restore
`.Random.seed` on exit - `.g_theme()` — Custom ggplot2 theme (serif,
white background) - `.permute_treatment()` — Handles cluster-level and
within-block permutation - `.validate_clusters_blocks()` /
`.validate_clusters_treatment()` — Input validation for experimental
designs

## Key Design Patterns

- **RNG discipline**: Every exported function saves `.Random.seed` on
  entry and restores it
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html). Internal seeds
  default to 1995.
- **Factor handling for grf**: Covariates go through a conversion
  pipeline — ordered factors become numeric, unordered factors get
  one-hot encoded via `model.matrix(~ . - 1, ...)`. This happens in both
  [`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md)
  and
  [`random_check()`](https://cetialphafive.github.io/MLbalance/reference/random_check.md).
- **Classifier backends in fastcpt**: Pluggable via `.gettrainmethod()`
  / `.gettestmethod()` factory functions. Each backend returns a train
  function and a predict function. `glmnet2` includes 2-way interactions
  and is binary-only.
- **Multi-arm treatment**: Joint K-class CPT on full data, then pairwise
  binary comparisons vs control for estimation. The `control` argument
  determines the reference level.
- **Overlap diagnostics**: When propensity scores are \< 0.05 or \>
  0.95, overlap-weighted estimates (Li, Morgan & Zaslavsky, 2018) are
  automatically computed via
  `grf::average_treatment_effect(target.sample = "overlap")`.

## Dependencies

**Imports**: distributional, estimatr, grf, ggdist, ggplot2, ranger,
rFerns **Suggests**: patchwork, mirai, glmnet, testthat

`grf` is the heaviest dependency — used for boosted regression forests
(propensity/outcome) and causal forests (ATE estimation).

## Code Conventions

- Assignment: `<-`, pipes: `|>`, variable names: `period.separated`
- Internal helpers are prefixed with `.` (e.g., `.g_theme`,
  `.permute_treatment`)
- roxygen2 with markdown enabled; NAMESPACE auto-generated
- `exportPattern("^[[:alpha:]]+")` is NOT used — exports are explicit in
  NAMESPACE
- testthat 3 edition
