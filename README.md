
# MLbalance <a href='https://github.com/samjfuller/MLbalance/blob/master/man/figures/mlbalance_sticker.png'><img src='man/figures/mlbalance_sticker.png' align="right" height="139" /></a>

MLbalance implements a novel machine learning balance test, the balance
permutation test, for experiments with binary, multiarm, and continuous
treatments. The purpose of this test is to detect failures of random
assignment, data fabrication, or imbalance across treatment arms. For more detail, see
[Rametta and Fuller (2024)](https://osf.io/preprints/osf/xcwt9).

This package is in development, any recommendations or comments welcome in the
issues section.

## Installation

Stable version is be available on Github here:

``` r
# install.packages("pak")
pak::pak("MLCause/MLbalance")
```

You can install the development version of MLbalance from
[GitHub](https://github.com/CetiAlphaFive/MLbalance) with:

``` r
# install.packages("pak")
pak::pak("CetiAlphaFive/MLbalance")
```

## Binary Treatment Example

Here is a basic example demonstrating the balance permutation test on a
simulated binary treatment DGP with multidimensional contamation of the
treatment assignment.

``` r
# install.packages("randomizr")
library(MLbalance)
#
set.seed(1995)
#
# Simple simulation 
n <- 1000
p <- 20
X <- matrix(rnorm(n*p,0,1),n,p)
w_real <- rbinom(n, 1, ifelse(.021 + abs(.4*X[,4] - .5*X[,8]) < 1, .021 + abs(.4*X[,4] - .5*X[,8]), 1))
# install.packages("randomizr")
w_sim <- randomizr::complete_ra(N = n,m = sum(w_real))
e <- rnorm(n,0,1)
y <- 2*w_real*X[,4] + 3*X[,2] -2*X[,8] + e
df <- data.frame(y,w_real,w_sim,X)
#
r.check <- random_check(W_real = df$w_real, #real treatment assignment vector 
                        W_sim  = df$w_sim, #simulated vector, comment out this argument to use permutated real assignment vector instead 
                        X      = subset(df,select = -c(y,w_real,w_sim)) #matrix of pretreatment covariates (or any covariates that SHOULD NOT be related to the assignment process/mechanism
             ); r.check$plot
#> Simulated Assignment Vector Provided, Null Distribution Generated Using Simulated Treatment Assignment.
#> 
#> 
#> Simple Count Table(s)
#> 
#> W_real
#>   0   1 
#> 520 480 
#> W_sim
#>   0   1 
#> 520 480 
#> 
#> 
#> Warning: Extreme Propensity scores detected (greater than .9 or less than .1).
#>                                       Examine $treat.props for more detail.
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
# to see variables important for predicting assignment, check r.check$imp.predictors 
```
