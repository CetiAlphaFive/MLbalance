# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

MLbalance (v0.2) provides ML-based covariate balance tests and causal effect estimation for experimental and observational data. The core tool is a fast classification permutation test (CPT) based on Gagnon-Bartsch & Shem-Tov (2019). If a classifier can distinguish treated from control units better than chance, balance fails. The package also estimates ATEs using four approaches: difference-in-means (DiM), inverse propensity weighted (IPW), outcome-adjusted, and AIPW (doubly robust) (all via `grf`).

## Repository Orientation

This dir (`package_exp/MLbalance/`) is the **active git repo and CRAN-track copy** — edit here. Do NOT confuse with sibling copies under the parent project: `package_push/MLbalance` (older, no test suite) and `package_push/MLbalance_exp/MLbalance`. The parent project root holds the paper + replication materials (separate `CLAUDE.md`), not package source.

## Build and Test Commands

```bash
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

Tests are slow (each `balance()` test fits multiple boosted forests). Use `perm.N = 50` and `n <= 250` in tests.

## Test Suite

testthat 3, six files under `tests/testthat/` (entry point `tests/testthat.R` → `test_check("MLbalance")`):
- `test-fastcpt.R` — engine: returned structure, p-value, null distribution.
- `test-fastcpt-backends.R` — optional backends (rpart/lda/qda), gated by Suggests + `requireNamespace`.
- `test-balance.R` — `balance()` with and without outcome `Y`; ATE estimators.
- `test-random_check.R` — `random_check()` structure and `vip()`.
- `test-coverage.R` — internal helpers (e.g. `.make_ci` confidence intervals).
- `test-adversarial.R` — hostile/edge inputs (e.g. character treatment labels, NAs in `X`, unequal group sizes).

## Architecture

Three-layer design, each layer an S3 class:

1. **`fastcpt(Z, T, ...)`** (`R/fastcpt.R`) — The permutation test engine. Trains classifiers (ferns, ranger forest, glmnet elastic net, linear probability, rpart, LDA, or QDA) on real treatment labels, then on `perm.N` permuted labels to build a null distribution. Returns p-value, test statistic, and null distribution. Supports `parallel = TRUE` via `mirai`. S3 methods in `R/fastcpt.plot.R`. The three new backends (rpart, lda, qda) are gated by Suggests + `requireNamespace`.

2. **`balance(Y, W, X, ...)`** (`R/balance.R`) — The main user-facing function. Calls `fastcpt()` for the balance test, then fits `grf::boosted_regression_forest` for propensity/outcome models and `grf::causal_forest` for ATE estimation (DiM, IPW, outcome-adjusted, AIPW). Handles multi-arm treatments via pairwise comparisons against a control level, with a joint K-class CPT. Detects extreme propensity scores and provides overlap-weighted (OW) fallback estimates. Returns S3 class "balance" with print/summary/plot methods.

3. **`random_check(W_real, X, ...)`** (`R/random_check.R`) — Lightweight diagnostic. Fits boosted RF propensity models on real vs. permuted/simulated treatment, returns overlapping propensity score distributions. Also exports `vip()` for variable importance from grf models.

**Shared utilities** (`R/utils.R`):
- `.save_rng_state()` / `.restore_rng_state()` — All exported functions save and restore `.Random.seed` on exit
- `.g_theme()` — Custom ggplot2 theme (serif, white background)
- `.permute_treatment()` — Handles cluster-level and within-block permutation
- `.validate_clusters_blocks()` / `.validate_clusters_treatment()` — Input validation for experimental designs

## Key Design Patterns

- **RNG discipline**: Every exported function saves `.Random.seed` on entry and restores it `on.exit()`. Internal seeds default to 1995.
- **Factor handling for grf**: Covariates go through a conversion pipeline — ordered factors become numeric, unordered factors get one-hot encoded via `model.matrix(~ . - 1, ...)`. This happens in both `balance()` and `random_check()`.
- **Classifier backends in fastcpt**: Pluggable via `.gettrainmethod()` / `.gettestmethod()` factory functions. Each backend returns a train function and a predict function. `glmnet2` includes 2-way interactions and is binary-only.
- **Multi-arm treatment**: Joint K-class CPT on full data, then pairwise binary comparisons vs control for estimation. The `control` argument determines the reference level.
- **Overlap diagnostics**: When propensity scores are < 0.05 or > 0.95, overlap-weighted estimates (Li, Morgan & Zaslavsky, 2018) are automatically computed via `grf::average_treatment_effect(target.sample = "overlap")`.
- **Propensity clamp**: pscores come from a *regression* forest, so they can exit [0,1] under strong covariate–treatment association. `balance()` clamps to `[eps, 1-eps]` (warning emitted, raw range reported) — without it IPW/AIPW go NaN. Separately, scores outside `overlap.threshold` (default `c(0.05, 0.95)`) trigger overlap-weighted (OW) fallback estimates. Return fields: `overlap_flag`, `n_extreme`, `overlap`.
- **Paired test**: `fastcpt(paired = TRUE)` errors on multi-class treatments and unequal group sizes (previously silently produced invalid permutations). `fastcpt3()` and the mlr3 backend were removed.

## Dependencies

**Imports**: estimatr, grf, ranger, rFerns
**Suggests**: distributional, ggdist, ggplot2, patchwork, mirai, glmnet, rpart, MASS, knitr, rmarkdown, tibble, withr, testthat

`grf` is the heaviest dependency — boosted regression forests (propensity/outcome) and causal forests (ATE). Plotting/distribution deps (ggplot2, ggdist, distributional, patchwork) are **Suggests, not Imports** — plot/summary code must guard them with `requireNamespace()`. Classifier backends (glmnet, rpart, MASS) and the `mirai` parallel backend are optional too.

## Code Conventions

- Assignment: `<-`, pipes: `|>`, variable names: `period.separated`
- Internal helpers are prefixed with `.` (e.g., `.g_theme`, `.permute_treatment`)
- roxygen2 with markdown enabled; NAMESPACE auto-generated
- `exportPattern("^[[:alpha:]]+")`  is NOT used — exports are explicit in NAMESPACE
- testthat 3 edition

## Docs & Release

- Vignette: `vignettes/MLbalance-workflow.Rmd` (knitr); long-form `balance.explainer.Rmd` at pkg root.
- pkgdown site: `_pkgdown.yml` → `docs/`. Rebuild with `pkgdown::build_site()`.
- `NEWS.md` tracks changes; `cran-comments.md` is the CRAN submission cover. Keep `R CMD check` CRAN-clean.
