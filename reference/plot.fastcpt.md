# Plot fastcpt Results

Creates a visualization of the Classification Permutation Test results.
Displays the permuted null distribution as a histogram with the observed
test statistic overlaid as a vertical line. The plot indicates whether
the test passed (p \> 0.05) or failed.

## Usage

``` r
# S3 method for class 'fastcpt'
plot(x, breaks = 25, ...)
```

## Arguments

- x:

  A fastcpt result object as returned by
  [`fastcpt`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md).

- breaks:

  Number of breaks for the histogram. Default is 25.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot object displaying the null distribution histogram with the
observed test statistic.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
n <- 200
p <- 10
Z <- matrix(rnorm(n * p), n, p)
T <- rep(c(1, 2), each = n/2)

# Run fast classification permutation test
result <- fastcpt(Z, T, class.methods = "forest", perm.N = 100)

# Plot the results
plot(result)
} # }
```
