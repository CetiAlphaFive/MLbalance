# Summary of fastcpt Results

Prints a formatted summary of the Classification Permutation Test
results including p-values, test statistics, and pass/fail status.

## Usage

``` r
# S3 method for class 'fastcpt'
summary(object, ...)
```

## Arguments

- object:

  A fastcpt result object as returned by
  [`fastcpt`](https://cetialphafive.github.io/MLbalance/reference/fastcpt.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the original object.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- fastcpt(Z, T, class.methods = "forest", perm.N = 100)
summary(result)
} # }
```
