# Suppress R CMD check notes for non-standard evaluation in internal functions
utils::globalVariables(c("pkgs"))

#' (fast) Classification Permutation Test
#'
#' @description
#' Code derived from the cpt package by Johann Gagnon-Bartsch. This version is optimized for speed, but the core structure is the same. The original package uses the RandomForest package, which is quite slow and lacks a variety of features available in other packages. This code converts to ranger and also adds an option for random ferns. Although obscure, random ferns are actually perfect for this - a very solid classifier that natively handles interactions and runs extremely fast on pretty much any hardware. I've also worked to add a version of logit that's faster using the RcppNumerical package, and added a basic parallel backend so the function runs reasonably quickly on larger datasets. Taken together these changes allow the user to run cpt in seconds rather than minutes or even hours for many datasets.
#'
#'Description of cpt: Non-parametric test for equality of multivariate distributions. Trains a classifier to classify (multivariate) observations as coming from one of several distributions. If the classifier is able to classify the observations better than would be expected by chance (using permutation inference), then the null hypothesis that the distributions are equal is rejected.
#'
#' @param Z The data. An n by p matrix, where n is the number of observations, and p is the number of covariates.
#' @param T The treatment variable. Is converted to a factor.
#' @param leaveout The number of observations from each treatment group to include in the test set. If 0, no data is left out and the in-sample test statistic is used. (See note below.) If an integer greater than or equal to 1, the number of observations from each treatment group to leave out. Values between 0 and 1 are converted to \code{ceiling(min(table(T))*leaveout)}.
#' @param class.methods A character vector of the different classification methods to use. Can be "forest", "ferns", or "logistic2fast". Default is "ferns" which is fast and handles interactions well.
#' @param metric Which test statistic to use. Can be "rate" (classification accuracy rate), "mse", or "probability". The default value ("probability") is recommended.
#' @param ensemble.metric Which test statistic to use for an ensemble classifier composed of all of the individual classifiers. Can be "vote", "mean.mse", or "mean.prob". The default value ("mean.prob") is recommended.
#' @param paired Do a paired permutation test. The data Z must be ordered such that the first observation with T==1 is paired with the first observation with T==2, the second observation with T==1 is paired with the second observation with T==2, etc. This can be accomplished by either letting the first n/2 rows be the treatment observations, and last n/2 rows being the control observations (in the same order), or by using the first two rows for the first pair, the second two rows for the second pair, etc.
#' @param perm.N The number of permutations.
#' @param leaveout.N The number of training set / test set iterations. In each iteration, a random test set is generated. Thus, test sets will typically overlap somewhat. There is one exception: If leaveout = 1 and leaveout.N = n, then a traditional leave-one-out procedure is used (each observation is left out exactly once).
#' @param comb.methods Which of the classifiers to include in the combined, overall p-value. Can be any subset of the classifiers specified in \code{class.methods} in addition to "ensemble" for the ensemble classifier.
#' @param comb.method The method for combining p-values from the individual classifiers in order to produce an overall p-value. The default ("fisher") is recommended. The other possible option is "min" which uses the minimum p-value. Note that in both cases, the combined p-value itself is not returned; rather, the combined p-value is treated as a test statistic, which is itself then subject to permutation analysis; what is returned is the resulting p-value from this analysis.
#' @param R.seed Random seed for R's set.seed (for reproducibility).
#' @param ranger.seed Random seed for ranger (for reproducibility).
#' @param parallel Logical. If TRUE, uses mirai for parallel processing across permutations.
#' @param alpha Significance level for pass/fail determination. Default is 0.05, which is appropriate for experimental contexts where even slight classification ability is concerning. For observational studies, a lower threshold may be more appropriate.
#' @param progress Logical. If TRUE (default), displays a progress bar during permutation testing.
#'
#' @return A list containing:
#' \item{pval}{The overall p-value, after combining results from the individual classifiers.}
#' \item{teststat}{The observed test statistics of the individual classifiers.}
#' \item{nulldist}{The permutation distributions of the individual classifiers.}
#' \item{pvals}{The p-values of the individual classifiers.}
#' \item{alpha}{The significance level used for pass/fail determination.}
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' n <- 200
#' p <- 10
#' Z <- matrix(rnorm(n * p), n, p)
#' T <- rep(c(1, 2), each = n/2)
#'
#' # Run fast classification permutation test
#' result <- fastcpt(Z, T, class.methods = "forest", perm.N = 100)
#' result$pval
#' }
#'
#' @export
#' @keywords multivariate
fastcpt <-
function (Z, T, leaveout = 0, class.methods = "ferns", metric = "probability",
    ensemble.metric = "mean.prob", paired = FALSE, perm.N = 1000,
    leaveout.N = 100, comb.methods = c(class.methods, "ensemble"),
    comb.method = "fisher", R.seed = 1995, ranger.seed = 1995, parallel = FALSE,
    alpha = 0.05, progress = TRUE)
{
    # Check required packages based on what's needed
    if (parallel && !requireNamespace("mirai", quietly = TRUE)) {
        stop("Package 'mirai' is required for parallel processing.", call. = FALSE)
    }
    if ("ferns" %in% class.methods && !requireNamespace("rFerns", quietly = TRUE)) {
        stop("Package 'rFerns' is required for ferns classifier.", call. = FALSE)
    }
    if ("forest" %in% class.methods) {
        if (!requireNamespace("ranger", quietly = TRUE)) {
            stop("Package 'ranger' is required for forest classifier.", call. = FALSE)
        }
    }
    if ("logistic2fast" %in% class.methods && !requireNamespace("RcppNumerical", quietly = TRUE)) {
        stop("Package 'RcppNumerical' is required for logistic2fast classifier.", call. = FALSE)
    }

    T = as.factor(T)
    train.methods = gettrainmethods(class.methods)
    test.methods = gettestmethods(class.methods)
    if (is.character(metric))
        metric = getmetric(metric)
    if (is.character(ensemble.metric))
        ensemble.metric = getensemble.metric(ensemble.metric)
    if (is.character(comb.method))
        comb.method = getcombmethod(comb.method)
    if ((leaveout > 0) && (leaveout < 1))
        leaveout = ceiling(min(table(T)) * leaveout)
    nulldist = matrix(NA, perm.N, length(class.methods) + 1)
    colnames(nulldist) = c(class.methods, "ensemble")

    # Get test statistic for real data
    teststat <- getteststat(Z, T, leaveout, train.methods, test.methods, metric, ensemble.metric, leaveout.N)
    
    if (paired) {
        T = as.numeric(T) - 1
        if (progress) pb <- utils::txtProgressBar(min = 0, max = perm.N, style = 3)
        for (i in 1:perm.N) {
            newT = T
            newT[T == 0] = stats::rbinom(length(T)/2, 1, 0.5)
            newT[T == 1] = 1 - newT[T == 0]
            newT = as.factor(newT)
            nulldist[i, ] = getteststat(Z, newT, leaveout, train.methods,
                test.methods, metric, ensemble.metric, leaveout.N)
            if (progress) utils::setTxtProgressBar(pb, i)
        }
        if (progress) close(pb)
    }
    else {
        # Pre-generate all permutations and per-iteration seeds (ensures reproducibility)
        base_seed <- sample.int(.Machine$integer.max, 1)
        perm_Ts <- lapply(1:perm.N, function(i) T[sample(length(T))])
        perm_seeds <- base_seed + (1:perm.N)

        if (parallel) {
            # Chunk permutations across workers (1 chunk per core, not 1 perm per task)
            n_workers <- parallel::detectCores() - 1
            chunk_ids <- split(1:perm.N, cut(1:perm.N, n_workers, labels = FALSE))

            mirai::daemons(n_workers, dispatcher = FALSE)
            on.exit(mirai::daemons(0), add = TRUE)

            # Each worker processes a chunk of pre-generated permutations
            results <- mirai::mirai_map(chunk_ids, function(ids) {
                # Load packages in worker
                for (pkg in pkgs) require(pkg, character.only = TRUE, quietly = TRUE)

                out <- matrix(NA, length(ids), length(train.methods) + 1)
                for (j in seq_along(ids)) {
                    set.seed(perm_seeds[ids[j]])
                    out[j, ] <- getteststat(Z, perm_Ts[[ids[j]]], leaveout, train.methods,
                        test.methods, metric, ensemble.metric, leaveout.N)
                }
                out
            }, Z = Z, leaveout = leaveout, train.methods = train.methods,
               test.methods = test.methods, metric = metric,
               ensemble.metric = ensemble.metric, leaveout.N = leaveout.N,
               getteststat = getteststat, train = train,
               applyclassifiers = applyclassifiers, softmax = softmax,
               pkgs = .packages(), perm_Ts = perm_Ts, perm_seeds = perm_seeds)[]

            # Reassemble results
            for (k in seq_along(chunk_ids)) {
                nulldist[chunk_ids[[k]], ] <- results[[k]]
            }
        } else {
            if (progress) pb <- utils::txtProgressBar(min = 0, max = perm.N, style = 3)
            for (i in 1:perm.N) {
                set.seed(perm_seeds[i])
                nulldist[i, ] = getteststat(Z, perm_Ts[[i]], leaveout, train.methods,
                    test.methods, metric, ensemble.metric, leaveout.N)
                if (progress) utils::setTxtProgressBar(pb, i)
            }
            if (progress) close(pb)
        }
    }
    pvals = rep(NA, ncol(nulldist))
    names(pvals) = colnames(nulldist)
    nullpvaldist = matrix(NA, perm.N, ncol(nulldist))
    colnames(nullpvaldist) = names(pvals)
    # Compute p-values with +1 correction to avoid p = 0 (Phipson & Smyth 2010)
    for (method.i in 1:ncol(nulldist)) {
        pvals[method.i] = (sum(nulldist[, method.i] >= teststat[method.i]) + 1) / (perm.N + 1)
        nullpvaldist[, method.i] = 1 - (rank(nulldist[, method.i],
            ties.method = "min") - 1)/perm.N
    }
    nullcombpvaldist = apply(nullpvaldist[, comb.methods, drop = FALSE],
        1, comb.method)
    pval = (sum(nullcombpvaldist <= comb.method(pvals[comb.methods])) + 1) / (perm.N + 1)

    if (length(class.methods) == 1) {
        result <- list(pval = pvals[1], teststat = teststat[1],
            nulldist = nulldist[, 1], pvals = pvals[1],
            class.methods = class.methods, metric = metric, perm.N = perm.N,
            alpha = alpha)
    } else {
        result <- list(pval = pval, teststat = teststat, nulldist = nulldist,
            pvals = pvals,
            class.methods = class.methods, metric = metric, perm.N = perm.N,
            alpha = alpha)
    }
    class(result) <- "fastcpt"
    return(result)
}

#' @keywords internal
#' @noRd
applyclassifiers <-
function (tstZ, tstT, classifiers, test.methods, metric, ensemble.metric,
    testistrain = FALSE)
{
    rval = rep(NA, length(classifiers) + 1)
    if (nrow(tstZ) == 1) {
        class.output = rep(NA, length(classifiers) * length(levels(tstT)))
        dim(class.output) = c(length(classifiers), length(levels(tstT)))
        for (i in 1:length(classifiers)) class.output[i, ] = test.methods[[i]](tstZ,
            classifiers[[i]], testistrain = testistrain)
        for (i in 1:length(classifiers)) rval[i] = metric(class.output[i,
            , drop = FALSE], tstT)
        dim(class.output) = c(1, dim(class.output))
        rval[length(classifiers) + 1] = ensemble.metric(class.output,
            tstT)
    }
    else {
        class.output = rep(NA, nrow(tstZ) * length(classifiers) *
            length(levels(tstT)))
        dim(class.output) = c(nrow(tstZ), length(classifiers),
            length(levels(tstT)))
        for (i in 1:length(classifiers)) class.output[, i, ] = test.methods[[i]](tstZ,
            classifiers[[i]], testistrain = testistrain)
        for (i in 1:length(classifiers)) rval[i] = metric(class.output[,
            i, ], tstT)
        rval[length(classifiers) + 1] = ensemble.metric(class.output,
            tstT)
    }
    return(rval)
}

#' @keywords internal
#' @noRd
getcombmethod <-
function (comb.method)
{
    if (comb.method == "fisher") {
        rval = function(x) {
            return(mean(log(x)))
        }
    }
    if (comb.method == "min") {
        rval = min
    }
    return(rval)
}

#' @keywords internal
#' @noRd
getensemble.metric <-
function (ensemble.metric)
{
    if (ensemble.metric == "vote") {
        rval = function(class.output, tstT) {
            temp = (class.output - as.vector(apply(class.output,
                c(1, 2), max))) == 0
            temp = temp/as.vector(apply(temp, c(1, 2), sum))
            votemat = apply(temp, c(1, 3), sum)
            indexmat = cbind(1:nrow(votemat), tstT)
            temp = votemat - apply(votemat, 1, max) == 0
            temp = temp/apply(temp, 1, sum)
            return(mean(votemat[indexmat]))
        }
    }
    if (ensemble.metric == "mean.prob") {
        rval = function(class.output, tstT) {
            meanprob = apply(class.output, c(1, 3), mean)
            indexmat = cbind(1:nrow(meanprob), tstT)
            return(mean(meanprob[indexmat]))
        }
    }
    if (ensemble.metric == "mean.log") {
        rval = function(class.output, tstT) {
            meanprob = apply(class.output, c(1, 3), mean)
            indexmat = cbind(1:nrow(meanprob), tstT)
            return(mean(log(meanprob[indexmat] + 1e-04)))
        }
    }
    if (ensemble.metric == "mean.mse") {
        rval = function(class.output, tstT) {
            meanprob = apply(class.output, c(1, 3), mean)
            indexmat = cbind(1:nrow(meanprob), tstT)
            meanprob[indexmat] = 1 - meanprob[indexmat]
            return(-mean(meanprob^2))
        }
    }
    return(rval)
}

#' @keywords internal
#' @noRd
getmetric <-
function (metric)
{
    if (metric == "probability") {
        rval = function(prob, tstT) {
            indexmat = cbind(1:nrow(prob), tstT)
            return(mean(prob[indexmat]))
        }
    }
    if (metric == "logscore") {
        rval = function(prob, tstT) {
            indexmat = cbind(1:nrow(prob), tstT)
            return(mean(log(prob[indexmat] + 1e-04)))
        }
    }
    else if (metric == "rate") {
        rval = function(prob, tstT) {
            indexmat = cbind(1:nrow(prob), tstT)
            temp = prob - apply(prob, 1, max) == 0
            temp = temp/apply(temp, 1, sum)
            return(mean(temp[indexmat]))
        }
    }
    else if (metric == "mse") {
        rval = function(prob, tstT) {
            indexmat = cbind(1:nrow(prob), tstT)
            prob[indexmat] = 1 - prob[indexmat]
            return(-mean(prob^2))
        }
    }
    return(rval)
}

#' @keywords internal
#' @noRd
gettestmethod <-
function (method)
{
    if (method == "forest") {
        rval = function(Z, classifier, testistrain = FALSE) {
            if (testistrain)
                return(classifier$predictions)
            else return(stats::predict(classifier, Z)$predictions)
        }
    }
    else if (method == "ferns") {
        rval = function(Z, classifier, testistrain = FALSE) {
            if (testistrain)
                return(softmax(classifier))
            else return(softmax(classifier, Z))
        }
    }
    else if (method == "logistic2fast") {
        rval = function(Z, classifier, testistrain = FALSE) {
            X <- stats::model.matrix(~.^2, data = as.data.frame(Z))

            # ensure test matrix matches training columns
            miss <- setdiff(classifier$cols, colnames(X))
            if (length(miss)) {
                X <- cbind(X, matrix(0, nrow(X), length(miss),
                                     dimnames = list(NULL, miss)))
            }
            X <- X[, classifier$cols, drop = FALSE]

            eta <- drop(X %*% classifier$beta)
            p1  <- 1 / (1 + exp(-eta))
            cbind(1 - p1, p1)  # n x 2 matrix, columns follow factor level order
        }
    }
    else {
        stop("Unknown classification method: ", method, call. = FALSE)
    }
    return(rval)
}

#' @keywords internal
#' @noRd
gettestmethods <-
function (class.methods)
{
    test.methods = list()
    for (i in 1:length(class.methods)) test.methods[[i]] = gettestmethod(class.methods[i])
    return(test.methods)
}

#' @keywords internal
#' @noRd
getteststat <-
function (Z, T, leaveout, train.methods, test.methods, metric,
    ensemble.metric, leaveout.N)
{
    if (leaveout == 0) {
        classifiers = train(Z, T, train.methods)
        return(applyclassifiers(Z, T, classifiers, test.methods,
            metric, ensemble.metric, testistrain = TRUE))
    }
    else if ((leaveout == 1) & (leaveout.N == nrow(Z))) {
        metric.mat = matrix(NA, leaveout.N, length(train.methods) +
            1)
        for (leaveout.i in 1:leaveout.N) {
            trnZ = Z[-leaveout.i, , drop = FALSE]
            trnT = T[-leaveout.i]
            tstZ = Z[leaveout.i, , drop = FALSE]
            tstT = T[leaveout.i]
            classifiers = train(trnZ, trnT, train.methods)
            metric.mat[leaveout.i, ] = applyclassifiers(tstZ,
                tstT, classifiers, test.methods, metric, ensemble.metric)
        }
        return(apply(metric.mat, 2, mean))
    }
    else {
        metric.mat = matrix(NA, leaveout.N, length(train.methods) +
            1)
        for (leaveout.i in 1:leaveout.N) {
            testset = rep(FALSE, length(T))
            for (i in 1:length(levels(T))) testset[sample(which(levels(T)[i] ==
                T), leaveout)] = TRUE
            trnZ = Z[!testset, , drop = FALSE]
            trnT = T[!testset]
            tstZ = Z[testset, , drop = FALSE]
            tstT = T[testset]
            classifiers = train(trnZ, trnT, train.methods)
            metric.mat[leaveout.i, ] = applyclassifiers(tstZ,
                tstT, classifiers, test.methods, metric, ensemble.metric)
        }
        return(apply(metric.mat, 2, mean))
    }
}

#' @keywords internal
#' @noRd
gettrainmethod <-
function (method)
{
    if (method == "forest") {
        rval = function(Z, T) {
            return(ranger::ranger(T ~ ., data = data.frame(T = T, Z), probability = TRUE, num.trees = 500L, num.threads = parallel::detectCores()-1))
        }
    }
    else if (method == "ferns") {
        rval = function(Z, T) {
            return(rFerns::rFerns(T ~ ., importance = "none", ferns = 1000L, depth = 5L, data = data.frame(T = T, Z)))
        }
    }
    else if (method == "logistic2fast") {
        rval = function(Z, T) {
            if (length(levels(T)) != 2)
                stop("logistic2fast supports binary T only.")
            X <- stats::model.matrix(~.^2, data = as.data.frame(Z))  # includes intercept
            y <- as.integer(T) - 1L                            # 0/1
            fit <- RcppNumerical::fastLR(x = X, y = y)
            beta <- if (!is.null(fit$coefficients)) as.numeric(fit$coefficients)
            else if (!is.null(fit$coef))     as.numeric(fit$coef)
            else stop("fastLR object lacks coefficients")
            list(beta = beta, cols = colnames(X))
        }
    }
    else {
        stop("Unknown classification method: ", method, call. = FALSE)
    }
    return(rval)
}

#' @keywords internal
#' @noRd
gettrainmethods <-
function (class.methods)
{
    train.methods = list()
    for (i in 1:length(class.methods)) train.methods[[i]] = gettrainmethod(class.methods[i])
    return(train.methods)
}

#' @keywords internal
#' @noRd
train <-
function (trnZ, trnT, train.methods)
{
    rval = list()
    for (i in 1:length(train.methods)) rval[[i]] = train.methods[[i]](trnZ,
        trnT)
    return(rval)
}

#' @keywords internal
#' @noRd
softmax <- function(fern.fit, newdata = NULL){
    if (is.null(newdata)) {
        S <- stats::predict(fern.fit, scores = TRUE)
    } else {
        S <- stats::predict(fern.fit, newdata, scores = TRUE)
    }
    S2 <- S - apply(S, 1, max)
    probs <- exp(S2)
    probs <- as.matrix(probs / rowSums(probs))
    return(probs)
}
