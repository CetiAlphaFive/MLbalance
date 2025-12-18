# Variable Importance Function

This convenience function takes a trained grf model object and returns a
data frame of variable importance scores using grf's simple variable
importance metric.

## Usage

``` r
vip(model)
```

## Arguments

- model:

  Trained GRF Model Object

## Examples

``` r
x <- data.frame(X1 = rnorm(1000))
y <- rnorm(1000)
model <- grf::regression_forest(X = x,Y = y)
vip(model)
#>   varname vip
#> 1      X1   1
```
