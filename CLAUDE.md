# MLbalance — CLAUDE.md

R package for machine learning balance tests and causal estimation in
experimental/observational data. Detects failures of random assignment, data
fabrication, or covariate imbalance. Authors: Jack Rametta, Sam Fuller.

## Quick Commands

```r
devtools::load_all()    # load package for interactive dev
devtools::document()    # regenerate NAMESPACE + man pages
devtools::test()        # run testthat suite
devtools::check()       # full R CMD check
pkgdown::build_site()   # rebuild docs site
devtools::build_readme() # reknit README.Rmd -> README.md
```

## Architecture

Three-layer design:

1. **`balance()`** (`R/balance.R`) — High-level wrapper. Runs balance test via
   `fastcpt()`, fits honest boosted RF propensity models via `grf`, computes
   four treatment effect estimators (DiM, IPW, outcome-adjusted, AIPW) with
   influence-function covariance. Supports binary and multi-arm treatments,
   cluster-randomized and blocked designs. Returns S3 class `"balance"`.

2. **`fastcpt()`** (`R/fastcpt.R`) — Standalone classification permutation
   test engine. Three classifier backends: random ferns (default), random
   forest, fast logistic. Cluster/block-aware permutation. Optional parallel
   via `mirai`. Returns S3 class `"fastcpt"`.

3. **`random_check()`** (`R/random_check.R`) — Diagnostic tool comparing real
   vs. null propensity score distributions from boosted RFs. Visual diagnostic
   without formal p-value. Also exports `vip()` for variable importance.

## File Layout

```
R/
├── MLbalance-package.R   # package-level docs (@keywords internal)
├── balance.R             # balance() + print/summary/plot S3 methods
├── fastcpt.R             # fastcpt() + internal classifier helpers
├── fastcpt.plot.R        # print/summary/plot S3 methods for fastcpt
├── random_check.R        # random_check() + vip()
└── utils.R               # shared internal helpers (.g_theme, .make_ci, etc.)
```

## Conventions

- **Internal functions**: dot-prefixed (`.compute_weights`), tagged `@keywords
  internal` and `@noRd`
- **Variable naming**: period.separated (`n.obs`, `avg.impurity`)
- **Assignment**: `<-` only
- **Pipes**: native `|>` (R 4.1+)
- **Default seed**: `1995` everywhere
- **RNG discipline**: all exported functions save/restore `.Random.seed` via
  `on.exit()` so they don't disturb the caller's RNG state
- **NSE globals**: declared with `utils::globalVariables()` at top of each file
- **roxygen2**: markdown enabled (`Roxygen: list(markdown = TRUE)`)
- **S3 method docs**: bundled under parent function via `@rdname`
- **Examples**: wrapped in `\donttest{}` (not `\dontrun{}`)

## Dependencies

**Imports:** distributional, estimatr, grf, ggdist, ggplot2, ranger, rFerns
**Suggests:** patchwork, mirai, RcppNumerical, testthat (>= 3.0.0)

Guard suggested packages with `requireNamespace()` checks at runtime.

## Testing

- testthat edition 3
- Tests in `tests/testthat/test-{balance,fastcpt,random_check}.R`
- All tests use small data (n=250, p=5), reduced `perm.N = 50`, `progress = FALSE`
- No `skip_on_cran()` — all tests run everywhere
- Run with: `devtools::test()` or `Rscript -e 'devtools::test()'`

## CI/CD

GitHub Actions in `.github/workflows/`:
- `check-standard.yaml` — R CMD check on macOS, Windows, Ubuntu (devel/release/oldrel)
- `lint.yaml` — `lintr::lint_package()`
- `pkgdown.yaml` — builds and deploys to `gh-pages`
- `test-coverage.yaml` — `covr` + Codecov

## Key Design Decisions

- `balance()` returns a rich S3 object with `$passed` for quick programmatic
  access and `plot()` for visual diagnostics
- Multi-arm treatments handled by pairwise comparisons against a `control` level
- Cluster-aware permutation: treatment must be constant within clusters;
  permutation happens at the cluster level
- Block-aware permutation: treatment is permuted within blocks
- Overlap diagnostics: when propensity scores are extreme (>0.975 or <0.025),
  overlap-weighted (OW) estimates are computed as fallback
- The shared ggplot theme (`.g_theme()` in utils.R) uses serif font for a
  publication-ready look
