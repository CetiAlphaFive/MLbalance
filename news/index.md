# Changelog

## MLbalance 0.2.1

- [`fastcpt()`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md)
  forest backend (`class.methods = "forest"`) now defaults to 100
  extremely-randomized trees (`splitrule = "extratrees"`), ~5x faster
  than the previous 500-tree gini default with equivalent size and power
  (validated over 6 core + 3 realistic DGPs at `perm.N = 1000`). It
  auto-falls back to `gini` when the data contain `NA` (extratrees
  cannot handle missing values).
- Any
  [`ranger::ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  argument can now be forwarded via `classifier.args` (e.g. `splitrule`,
  `num.random.splits`, `min.node.size`, `sample.fraction`).
- `write.forest` is skipped in the out-of-bag path (`leaveout = 0`) for
  speed; the leave-out predict path now restores trained column names
  before prediction (fixes a pre-existing failure with unnamed covariate
  matrices).

## MLbalance 0.2

- Initial CRAN submission.
- [`balance()`](https://cetialphafive.github.io/MLbalance/reference/balance.md):
  unified function for covariate balance assessment and treatment effect
  estimation supporting binary and multi-arm treatments.
- [`fastcpt()`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md):
  fast classification permutation test with support for random ferns,
  ranger forests, and fast logistic regression.
- [`random_check()`](https://cetialphafive.github.io/MLbalance/reference/random_check.md):
  balance permutation test using boosted regression forests with
  diagnostic plots.
- [`vip()`](https://cetialphafive.github.io/MLbalance/reference/vip.md):
  variable importance convenience function for grf model objects.
