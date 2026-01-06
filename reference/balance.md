# Complete Balance Assessment and Treatment Effect Estimation

A unified function for covariate balance assessment and treatment effect
estimation. Combines a formal balance test (via classification
permutation test), visual diagnostics (propensity score distributions),
and treatment effect estimates using both difference-in-means and
augmented inverse propensity weighting (AIPW).

Supports both binary and multi-arm treatments. For multi-arm treatments,
pairwise comparisons are made between each treatment arm and the control
group.

## Usage

``` r
balance(
  Y = NULL,
  W,
  X,
  alpha = 0.05,
  perm.N = 1000,
  class.method = "ferns",
  seed = 1995,
  control = NULL,
  fastcpt.args = list()
)

# S3 method for class 'balance'
print(x, ...)

# S3 method for class 'balance'
summary(object, ...)

# S3 method for class 'balance'
plot(x, which = "all", combined = TRUE, breaks = 15, ...)
```

## Arguments

- Y:

  Outcome vector (numeric) or `NULL`. If `NULL`, treatment effect
  estimation (and the treatment effect plot) is skipped.

- W:

  Treatment assignment vector. Can be binary (0/1, logical) or multi-arm
  (factor, character, or integer with \>2 levels).

- X:

  Pre-treatment covariate matrix or data frame.

- alpha:

  Significance level for balance test. Default is 0.05.

- perm.N:

  Number of permutations for the balance test. Default is 1000.

- class.method:

  Classification method for balance test. Default is "ferns".

- seed:

  Random seed for reproducibility. Default is 1995.

- control:

  Optional. The value in `W` to use as the control group. If `NULL`
  (default), the first factor level is used as control. A message is
  displayed indicating the control assumption.

- fastcpt.args:

  A named list of additional arguments to pass to
  [`fastcpt`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md).
  For example, `fastcpt.args = list(parallel = TRUE, leaveout = 0.2)`.

- x:

  A balance result object.

- ...:

  Additional arguments (currently unused).

- object:

  A balance result object (for summary method).

- which:

  Character vector specifying which plots to create. Options are
  "pscores", "null_dist", "effects", or "all".

- combined:

  Logical. If TRUE, displays all three plots in a combined panel.
  Default is TRUE.

- breaks:

  Number of breaks for histograms. Default is 15.

## Value

A list of class "balance" containing:

- balance_test:

  Results from fastcpt including p-value and propensity scores. For
  multi-arm, a named list with one entry per treatment arm.

- dim:

  Difference-in-means estimate with standard error and confidence
  interval (only if `Y` is provided). For multi-arm, a named list.

- aipw:

  Doubly robust estimate from causal forest (with propensity weighting)
  with standard error and CI (only if `Y` is provided). For multi-arm, a
  named list.

- aipw_const:

  Outcome-adjusted estimate from causal forest (no propensity weighting)
  with SE and CI (only if `Y` is provided). For multi-arm, a named list.

- passed:

  Logical indicating whether the balance test passed. For multi-arm, a
  named logical vector.

- alpha:

  The significance level used.

- cf:

  The fitted causal_forest object(s) for advanced users (only if `Y` is
  provided). For multi-arm, a named list.

- control:

  The control level used.

- arms:

  Character vector of treatment arm names (excluding control).

- multiarm:

  Logical indicating whether this is a multi-arm analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data (binary treatment)
n <- 500
p <- 10
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
Y <- W * 0.5 + X[,1] * 0.3 + rnorm(n)

# Run complete balance assessment
result <- balance(Y, W, X)
result
summary(result)
plot(result)

# Multi-arm example
W_multi <- sample(c("Control", "Treatment A", "Treatment B"), n, replace = TRUE)
Y_multi <- (W_multi == "Treatment A") * 0.3 + (W_multi == "Treatment B") * 0.6 + rnorm(n)
result_multi <- balance(Y_multi, W_multi, X, control = "Control")
plot(result_multi)
} # }
```
