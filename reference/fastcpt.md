# (fast) Classification Permutation Test

Code derived from the cpt package by Johann Gagnon-Bartsch. This version
is optimized for speed, but the core structure is the same. The original
package uses the RandomForest package, which is quite slow and lacks a
variety of features available in other packages. This code converts to
ranger and also adds an option for random ferns. Although obscure,
random ferns are great for this - a very solid classifier that natively
handles interactions and runs extremely fast on pretty much any
hardware. I've also worked to add a version of logit using glmnet with
elastic net regularization to prevent coefficient blowup under
separation, and added a basic parallel backend so the function runs
reasonably quickly on larger datasets. Taken together these changes
allow the user to run cpt in seconds rather than minutes or even hours
for many datasets.

Description of cpt: Non-parametric test for equality of multivariate
distributions. Trains a classifier to classify (multivariate)
observations as coming from one of several distributions. If the
classifier is able to classify the observations better than would be
expected by chance (using permutation inference), then the null
hypothesis that the distributions are equal is rejected.

## Usage

``` r
fastcpt(
  Z,
  T,
  leaveout = 0,
  class.methods = "ferns",
  metric = "probability",
  ensemble.metric = "mean.prob",
  paired = FALSE,
  clusters = NULL,
  blocks = NULL,
  perm.N = 1000,
  leaveout.N = 100,
  comb.methods = c(class.methods, "ensemble"),
  comb.method = "fisher",
  R.seed = 1995,
  ranger.seed = 1995,
  parallel = FALSE,
  alpha = 0.05,
  progress = interactive(),
  classifier.args = list()
)

# S3 method for class 'fastcpt'
plot(x, breaks = 25, ...)

# S3 method for class 'fastcpt'
summary(object, ...)

# S3 method for class 'fastcpt'
print(x, ...)
```

## Arguments

- Z:

  The data. An n by p matrix, where n is the number of observations, and
  p is the number of covariates.

- T:

  The treatment variable. Is converted to a factor.

- leaveout:

  The number of observations from each treatment group to include in the
  test set. If 0, no data is left out and the in-sample test statistic
  is used. (See note below.) If an integer greater than or equal to 1,
  the number of observations from each treatment group to leave out.
  Values between 0 and 1 are converted to
  `ceiling(min(table(T))*leaveout)`.

- class.methods:

  A character vector of the different classification methods to use. Can
  be "forest", "ferns", "glmnet2", or "lm". Default is "ferns" which is
  fast and handles interactions well.

- metric:

  Which test statistic to use. Can be "rate", "mse", "logscore", or
  "probability" (default, recommended).

- ensemble.metric:

  Which test statistic to use for an ensemble classifier composed of all
  of the individual classifiers. Can be "vote", "mean.mse", "mean.log",
  or "mean.prob" (default, recommended).

- paired:

  Do a paired permutation test. The data Z must be ordered such that the
  first observation with T==1 is paired with the first observation with
  T==2, the second observation with T==1 is paired with the second
  observation with T==2, etc. This can be accomplished by either letting
  the first n/2 rows be the treatment observations, and last n/2 rows
  being the control observations (in the same order), or by using the
  first two rows for the first pair, the second two rows for the second
  pair, etc.

- clusters:

  Optional vector of cluster identifiers (same length as `T`). When
  provided, permutations shuffle treatment labels at the cluster level
  rather than the individual level. Treatment must be constant within
  each cluster.

- blocks:

  Optional vector of block identifiers (same length as `T`). When
  provided, permutations are restricted to within each block. Cannot be
  used together with `paired`.

- perm.N:

  The number of permutations.

- leaveout.N:

  The number of training set / test set iterations. In each iteration, a
  random test set is generated. Thus, test sets will typically overlap
  somewhat. There is one exception: If leaveout = 1 and leaveout.N = n,
  then a traditional leave-one-out procedure is used (each observation
  is left out exactly once).

- comb.methods:

  Which of the classifiers to include in the combined, overall p-value.
  Can be any subset of the classifiers specified in `class.methods` in
  addition to "ensemble" for the ensemble classifier.

- comb.method:

  The method for combining p-values from the individual classifiers in
  order to produce an overall p-value. The default ("fisher") is
  recommended. The other possible option is "min" which uses the minimum
  p-value. Note that in both cases, the combined p-value itself is not
  returned; rather, the combined p-value is treated as a test statistic,
  which is itself then subject to permutation analysis; what is returned
  is the resulting p-value from this analysis. The "fisher" option
  computes `mean(log(p))`, a monotone transformation of the classical
  Fisher statistic; since comparison is against a permutation null
  rather than chi-squared, the scaling is immaterial.

- R.seed:

  Random seed for R's set.seed (for reproducibility).

- ranger.seed:

  Random seed for ranger (for reproducibility).

- parallel:

  Logical. If TRUE, uses mirai for parallel processing across
  permutations.

- alpha:

  Significance level for pass/fail determination. Default is 0.05, which
  is appropriate for experimental contexts where even slight
  classification ability is concerning. For observational studies, a
  lower threshold may be more appropriate.

- progress:

  Logical. If TRUE, displays a progress bar during permutation testing.
  Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- classifier.args:

  Optional named list of hyperparameters forwarded to the classifier
  training functions. Supported keys: `num.trees` (for ranger forest,
  default 500), `ferns` (number of ferns for rFerns, default 500),
  `depth` (fern depth for rFerns, default 5), `nfolds` (for cv.glmnet,
  default 5), `alpha` (elastic net mixing for cv.glmnet, default 0.5).

- x:

  A fastcpt result object (for plot and print methods).

- breaks:

  Number of breaks for the histogram. Default is 25.

- ...:

  Additional arguments (currently unused).

- object:

  A fastcpt result object (for summary method).

## Value

A list containing:

- pval:

  The overall p-value, after combining results from the individual
  classifiers.

- teststat:

  The observed test statistics of the individual classifiers.

- nulldist:

  The permutation distributions of the individual classifiers.

- pvals:

  The p-values of the individual classifiers.

- alpha:

  The significance level used for pass/fail determination.

- class.methods:

  Character vector of classification methods used.

- metric:

  The metric function used.

- metric_name:

  Character name of the metric.

- perm.N:

  Number of permutations.

- clusters:

  Cluster identifiers (if provided, otherwise `NULL`).

- blocks:

  Block identifiers (if provided, otherwise `NULL`).

## Examples

``` r
# \donttest{
# Generate example data
n <- 200
p <- 10
Z <- matrix(rnorm(n * p), n, p)
T <- rep(c(1, 2), each = n/2)

# Run fast classification permutation test
result <- fastcpt(Z, T, class.methods = "forest", perm.N = 100)
result$pval
#>    forest 
#> 0.2970297 
# }
```
