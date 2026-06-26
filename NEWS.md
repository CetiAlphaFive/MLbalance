# MLbalance (Unreleased)

* Removed `fastcpt3()` and the `mlr3` / `mlr3learners` dependency. Added native `rpart`, `lda`, and `qda` backends to `fastcpt()` (all optional Suggests).
* Bug fix in `fastcpt(paired = TRUE)`: now errors clearly on multi-class treatments and on unequal group sizes. Previously these cases silently produced invalid permutations (multi-class rows untouched; binary unequal-size groups recycled from `rbinom(length(T)/2, ...)`).
* `plot.balance()`: the propensity-score panel (Panel A) now compares propensity score distributions between treatment arms (control vs treated) instead of real vs permuted-null. The classification permutation test panel (B) and treatment-effect panel (C) are unchanged.

# MLbalance 0.2

* Initial CRAN submission.
* `balance()`: unified function for covariate balance assessment and treatment effect estimation supporting binary and multi-arm treatments.
* `fastcpt()`: fast classification permutation test with support for random ferns, ranger forests, and fast logistic regression.
* `random_check()`: balance permutation test using boosted regression forests with diagnostic plots.
* `vip()`: variable importance convenience function for grf model objects.
