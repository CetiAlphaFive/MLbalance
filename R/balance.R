# Suppress R CMD check notes for non-standard evaluation
utils::globalVariables(c("estimator", "estimate", "var", "val", ".dist", "arm", "comparison",
                         "mean_real", "mean_null", "testval", "pval", "status_text", "color",
                         "weight"))

#' Complete Balance Assessment and Treatment Effect Estimation
#'
#' @description
#' A unified function for covariate balance assessment and treatment effect estimation.
#' Combines a formal balance test (via classification permutation test), visual diagnostics
#' (propensity score distributions), and treatment effect estimates using both difference-in-means
#' and AIPW (augmented inverse propensity weighting).
#'
#' Supports both binary and multi-arm treatments. For multi-arm treatments, pairwise comparisons
#' are made between each treatment arm and the control group.
#'
#' @param Y Outcome vector (numeric) or \code{NULL}. If \code{NULL}, treatment effect
#'   estimation (and the treatment effect plot) is skipped.
#' @param W Treatment assignment vector. Can be binary (0/1, logical) or multi-arm (factor, character, or integer with >2 levels).
#' @param X Pre-treatment covariate matrix or data frame.
#' @param alpha Significance level for balance test. Default is 0.05.
#' @param perm.N Number of permutations for the balance test. Default is 1000.
#' @param class.method Classification method for balance test. Can be "ferns" (default),
#'   "forest", or "glmnet2". To use an ensemble of classifiers, pass
#'   \code{fastcpt.args = list(class.methods = c("ferns", "forest"))}.
#' @param seed Random seed for reproducibility. Default is 1995.
#' @param control Optional. The value in \code{W} to use as the control group. If \code{NULL} (default),
#'   the first factor level is used as control. A message is displayed indicating the control assumption.
#' @param clusters Optional vector of cluster identifiers (same length as \code{W}). When provided,
#'   permutations in the balance test shuffle treatment labels at the cluster level, and treatment
#'   effect standard errors use cluster-robust variance estimators. Treatment must be constant
#'   within each cluster.
#' @param blocks Optional vector of block identifiers (same length as \code{W}). When provided,
#'   permutations in the balance test are restricted to within each block.
#' @param num.trees Number of trees used in \code{grf::causal_forest()} for treatment
#'   effect estimation. Default is 2000.
#' @param overlap.threshold Numeric vector of length 2 giving the lower and upper propensity
#'   score thresholds for flagging overlap issues. Default is \code{c(0.05, 0.95)}.
#'   When any propensity scores fall outside these bounds, overlap-weighted estimates
#'   are automatically computed.
#' @param fastcpt.args A named list of additional arguments to pass to \code{\link{fastcpt}}.
#'   For example, \code{fastcpt.args = list(parallel = TRUE, leaveout = 0.2)}.
#'   You can also pass classifier-specific hyperparameters through this list, e.g.,
#'   \code{fastcpt.args = list(classifier.args = list(num.trees = 1000))} for ranger,
#'   \code{list(classifier.args = list(ferns = 1000))} for rFerns, or
#'   \code{list(classifier.args = list(nfolds = 10))} for cv.glmnet.
#'   You can also use this to run an ensemble of classifiers:
#'   \code{fastcpt.args = list(class.methods = c("ferns", "forest"))}.
#'
#' @return A list of class "balance" containing:
#' \item{balance_test}{Results from fastcpt including p-value and propensity scores. For multi-arm, a named list with one entry per treatment arm.}
#' \item{dim}{Difference-in-means estimate with standard error and confidence interval (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{ipw}{IPW estimate using propensity scores from the boosted regression forest, with SE and CI (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{aipw}{AIPW (doubly robust) estimate from causal forest (with propensity weighting) with standard error and CI (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{oadj}{Outcome-adjusted estimate from causal forest (no propensity weighting) with SE and CI (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{passed}{Logical indicating whether the balance test passed. For multi-arm, a named logical vector.}
#' \item{alpha}{The significance level used.}
#' \item{cf}{The fitted causal_forest object(s) for advanced users (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{imp.predictors}{Variable importance scores from the propensity model, computed via \code{\link{vip}}. For multi-arm, a named list.}
#' \item{control}{The control level used.}
#' \item{arms}{Character vector of treatment arm names (excluding control).}
#' \item{multiarm}{Logical indicating whether this is a multi-arm analysis.}
#' \item{overlap_flag}{Logical indicating whether overlap issues were detected.}
#' \item{overlap}{Overlap-weighted estimates (if \code{overlap_flag} is \code{TRUE}).}
#' \item{n_extreme}{Number of observations with extreme propensity scores.}
#' \item{pscores_real}{Propensity scores from the real treatment assignment.}
#' \item{pscores_null}{Propensity scores from a permuted treatment assignment.}
#' \item{n}{Number of observations.}
#' \item{n_treated}{Number of treated units (binary case).}
#' \item{n_control}{Number of control units (binary case).}
#' \item{n_per_arm}{Named vector of sample sizes per arm (multi-arm case).}
#' \item{clusters}{Cluster identifiers (if provided).}
#' \item{blocks}{Block identifiers (if provided).}
#' \item{ate_cov}{Covariance matrix of ATE estimates (used in divergence tests).}
#'
#' @examples
#' \donttest{
#' # Generate example data (binary treatment)
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- W * 0.5 + X[,1] * 0.3 + rnorm(n)
#'
#' # Run complete balance assessment
#' result <- balance(Y, W, X)
#' result
#' summary(result)
#' plot(result)
#'
#' # Multi-arm example
#' W_multi <- sample(c("Control", "Treatment A", "Treatment B"), n, replace = TRUE)
#' Y_multi <- (W_multi == "Treatment A") * 0.3 + (W_multi == "Treatment B") * 0.6 + rnorm(n)
#' result_multi <- balance(Y_multi, W_multi, X, control = "Control")
#' plot(result_multi)
#' }
#'
#' @export
balance <- function(Y = NULL, W, X, alpha = 0.05, perm.N = 1000, class.method = "ferns", seed = 1995,
                    control = NULL, clusters = NULL, blocks = NULL,
                    num.trees = 2000, overlap.threshold = c(0.05, 0.95),
                    fastcpt.args = list()) {

  # Save and restore RNG state
  old_seed <- .save_rng_state()
  on.exit(.restore_rng_state(old_seed), add = TRUE)

  # Check for grf
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required for propensity score modeling (and optional treatment effect estimation).", call. = FALSE)
  }

  # Input validation
  W_levels <- unique(W)
  n_arms <- length(W_levels)

  if (n_arms < 2) stop("W must have at least 2 unique values.", call. = FALSE)
  if (length(W) != nrow(X)) {
    stop("W and X must have the same number of observations.", call. = FALSE)
  }
  # NA/NaN checks
  if (any(is.na(W)))
    stop("W contains NA values.", call. = FALSE)
  if (any(is.na(as.matrix(X))))
    stop("X contains NA or NaN values. Remove or impute before running balance().", call. = FALSE)

  # Inf check — catches log(0), 1/0, etc.
  X_df_tmp <- as.data.frame(X)
  num_cols <- vapply(X_df_tmp, is.numeric, logical(1))
  if (any(num_cols) && any(is.infinite(as.matrix(X_df_tmp[num_cols]))))
    stop("X contains Inf or -Inf values. Remove or replace before running balance().", call. = FALSE)

  if (!is.null(Y)) {
    if (!is.numeric(Y)) stop("Y must be a numeric vector.", call. = FALSE)
    if (length(Y) != length(W)) {
      stop("Y and W must have the same number of observations.", call. = FALSE)
    }
    if (any(is.na(Y)))
      stop("Y contains NA values. Remove or impute before running balance().", call. = FALSE)
    if (any(is.infinite(Y)))
      stop("Y contains Inf or -Inf values. Remove or replace before running balance().", call. = FALSE)
  }

  .validate_clusters_blocks(clusters, blocks, length(W), paired = FALSE)
  if (!is.null(clusters)) .validate_clusters_treatment(W, clusters)

  # Determine control level
  W_factor <- as.factor(W)
  all_levels <- levels(W_factor)

  if (is.null(control)) {
    control <- all_levels[1]
    if (n_arms > 2 && getOption("MLbalance.verbose", TRUE)) {
      message(sprintf("Multi-arm treatment detected (%d arms). Using '%s' as control (first factor level).\n  To specify a different control, use the 'control' argument.",
                      n_arms, control))
    }
  } else {
    if (!(control %in% all_levels)) {
      stop(sprintf("Specified control '%s' not found in W. Available levels: %s",
                   control, paste(all_levels, collapse = ", ")), call. = FALSE)
    }
    # Reorder factor levels so control is first
    W_factor <- factor(W, levels = c(control, setdiff(all_levels, control)))
    all_levels <- levels(W_factor)
    if (getOption("MLbalance.verbose", TRUE)) {
      message(sprintf("Using '%s' as control group.", control))
    }
  }

  # Treatment arms (excluding control)
  treatment_arms <- setdiff(all_levels, control)
  multiarm <- n_arms > 2

  # Convert X to data frame; X_df is passed to fastcpt (ranger/ferns handle factors natively)
  X_df <- as.data.frame(X)

  # Convert Date/POSIXt columns to numeric (days since epoch) — grf and ferns need numeric/factor
  for (j in seq_along(X_df)) {
    if (inherits(X_df[[j]], "Date") || inherits(X_df[[j]], "POSIXt"))
      X_df[[j]] <- as.numeric(X_df[[j]])
  }

  # Build numeric matrix for grf: ordered factors → numeric, unordered → one-hot
  X_matrix <- .prepare_covariates_for_grf(X_df)
  if (ncol(X_matrix) > ncol(X_df) * 5) {
    warning("Factor expansion created many columns. Consider collapsing high-cardinality factors.",
            call. = FALSE)
  }

  # Set seed
  set.seed(seed)

  # Initialize storage for results
  balance_test_joint <- NULL
  balance_tests <- list()
  pscores_real_list <- list()
  pscores_null_list <- list()
  dim_results <- list()
  ipw_results <- list()
  aipw_results <- list()
  oadj_results <- list()
  cf_list <- list()
  ate_cov_list <- list()
  overlap_flag_list <- list()
  overlap_results <- list()
  n_extreme_list <- list()
  passed_vec <- c()
  n_per_arm <- list()
  vip_list <- list()

  # For multi-arm: single joint K-class balance test on full data
  if (multiarm) {
    if (class.method == "glmnet2") {
      stop("glmnet2 only supports binary treatment. Use 'ferns' or 'forest' for multi-arm.", call. = FALSE)
    }
    fastcpt_base_args <- list(
      Z = X_df, T = W_factor,
      class.methods = class.method,
      perm.N = perm.N, alpha = alpha, R.seed = seed,
      clusters = clusters, blocks = blocks
    )
    fastcpt_all_args <- utils::modifyList(fastcpt_base_args, fastcpt.args)
    balance_test_joint <- do.call(fastcpt, fastcpt_all_args)
    passed_joint <- balance_test_joint$pval > alpha
  }

  # Loop over each treatment arm (pairwise vs control)
  for (arm in treatment_arms) {
    arm_name <- if (multiarm) arm else "treated"

    # Subset to control + this arm
    idx <- W_factor %in% c(control, arm)
    W_binary <- as.integer(W_factor[idx] == arm)  # 1 = treatment arm, 0 = control
    X_sub <- X_df[idx, , drop = FALSE]
    X_matrix_sub <- X_matrix[idx, , drop = FALSE]
    Y_sub <- if (!is.null(Y)) Y[idx] else NULL
    clusters_sub <- if (!is.null(clusters)) clusters[idx] else NULL
    blocks_sub <- if (!is.null(blocks)) blocks[idx] else NULL

    # Augment covariate matrix with block dummies for grf (which has no blocks arg)
    if (!is.null(blocks_sub)) {
      block_dummies <- stats::model.matrix(~ as.factor(blocks_sub) - 1)
      if (ncol(block_dummies) > 1) {
        block_dummies <- block_dummies[, -1, drop = FALSE]
      } else {
        block_dummies <- NULL
      }
      if (!is.null(block_dummies)) {
        colnames(block_dummies) <- paste0(".block_", seq_len(ncol(block_dummies)))
        X_matrix_sub <- cbind(X_matrix_sub, block_dummies)
      }
    }

    n_per_arm[[arm_name]] <- list(
      control = sum(W_binary == 0),
      treated = sum(W_binary == 1)
    )

    # 1. Run fastcpt for formal balance test (binary only; multi-arm uses joint test above)
    if (!multiarm) {
      fastcpt_base_args <- list(
        Z = X_sub,
        T = W_binary,
        class.methods = class.method,
        perm.N = perm.N,
        alpha = alpha,
        R.seed = seed,
        clusters = clusters_sub,
        blocks = blocks_sub
      )
      fastcpt_all_args <- utils::modifyList(fastcpt_base_args, fastcpt.args)
      balance_test <- do.call(fastcpt, fastcpt_all_args)
      balance_tests[[arm_name]] <- balance_test
      passed_vec[arm_name] <- balance_test$pval > alpha
    }

    # 2. Fit honest propensity models for visualization (boosted RF from grf)
    n_sub <- length(W_binary)
    n_min_group <- min(sum(W_binary == 0), sum(W_binary == 1))

    .wrap_grf <- function(expr) {
      tryCatch(expr, error = function(e) {
        if (grepl("honesty fraction", e$message, fixed = TRUE)) {
          stop(sprintf(
            "Sample too small or too imbalanced for grf propensity model (n = %d, min group = %d in '%s' vs '%s'). Use fastcpt() directly for balance testing.",
            n_sub, n_min_group, arm, control), call. = FALSE)
        }
        stop(e)
      })
    }

    prop_real <- .wrap_grf(grf::boosted_regression_forest(
      X = X_matrix_sub,
      Y = W_binary,
      honesty = TRUE,
      tune.parameters = "all",
      seed = seed,
      clusters = clusters_sub
    ))
    pscores_real_list[[arm_name]] <- as.numeric(prop_real$predictions)

    # Variable importance from the propensity model
    vip_list[[arm_name]] <- vip(prop_real$forests[[1]])

    # Null (permuted) treatment propensities
    set.seed(seed + which(treatment_arms == arm))
    W_null <- .permute_treatment(W_binary, clusters_sub, blocks_sub)
    prop_null <- .wrap_grf(grf::boosted_regression_forest(
      X = X_matrix_sub,
      Y = W_null,
      honesty = TRUE,
      tune.parameters = "all",
      seed = seed,
      clusters = clusters_sub
    ))
    pscores_null_list[[arm_name]] <- as.numeric(prop_null$predictions)

    # Detect extreme propensity scores
    ps_real <- pscores_real_list[[arm_name]]
    n_extreme <- sum(ps_real < overlap.threshold[1]) + sum(ps_real > overlap.threshold[2])
    has_overlap_issue <- n_extreme > 0
    overlap_flag_list[[arm_name]] <- has_overlap_issue
    n_extreme_list[[arm_name]] <- n_extreme

    # 3. Optional treatment effect estimation (only if Y provided)
    if (!is.null(Y_sub)) {
      # Compute difference-in-means via estimatr (handles clusters automatically)
      dim_fit <- estimatr::difference_in_means(
        Y_sub ~ W_binary,
        clusters = clusters_sub,
        blocks = blocks_sub
      )
      dim_est <- dim_fit$coefficients[[1]]
      dim_se  <- dim_fit$std.error[[1]]

      dim_results[[arm_name]] <- list(
        estimate = dim_est,
        std.err = dim_se,
        ci = .make_ci(dim_est, dim_se)
      )

      # Needed for influence function below
      Y1 <- Y_sub[W_binary == 1]
      Y0 <- Y_sub[W_binary == 0]

      e_hat <- pscores_real_list[[arm_name]]

      # Shared outcome model (boosted RF), used by all outcome-adjusting forests
      outcome_forest <- .wrap_grf(grf::boosted_regression_forest(
        X = X_matrix_sub,
        Y = Y_sub,
        honesty = TRUE,
        tune.parameters = "all",
        seed = seed,
        clusters = clusters_sub
      ))
      Y_hat <- as.numeric(outcome_forest$predictions)

      # IPW: boosted-RF propensity, flat outcome model
      cf_ipw <- grf::causal_forest(
        X     = X_matrix_sub,
        Y     = Y_sub,
        W     = W_binary,
        W.hat = e_hat,
        Y.hat = rep(mean(Y_sub), length(Y_sub)),
        seed  = seed,
        num.trees = num.trees,
        clusters = clusters_sub
      )
      ate_ipw <- grf::average_treatment_effect(cf_ipw, target.sample = "all")
      ipw_results[[arm_name]] <- list(
        estimate = ate_ipw["estimate"],
        std.err  = ate_ipw["std.err"],
        ci       = .make_ci(ate_ipw["estimate"], ate_ipw["std.err"])
      )

      # Outcome-adjusted: constant propensity, boosted-RF outcome model
      cf_const <- grf::causal_forest(
        X     = X_matrix_sub,
        Y     = Y_sub,
        W     = W_binary,
        W.hat = rep(mean(W_binary), length(W_binary)),
        Y.hat = Y_hat,
        seed  = seed,
        num.trees = num.trees,
        clusters = clusters_sub
      )
      ate_const <- grf::average_treatment_effect(cf_const, target.sample = "all")
      oadj_results[[arm_name]] <- list(
        estimate = ate_const["estimate"],
        std.err  = ate_const["std.err"],
        ci       = .make_ci(ate_const["estimate"], ate_const["std.err"])
      )

      # AIPW: boosted-RF propensity + boosted-RF outcome model
      cf <- grf::causal_forest(
        X     = X_matrix_sub,
        Y     = Y_sub,
        W     = W_binary,
        W.hat = e_hat,
        Y.hat = Y_hat,
        seed  = seed,
        num.trees = num.trees,
        clusters = clusters_sub
      )
      cf_list[[arm_name]] <- cf
      ate <- grf::average_treatment_effect(cf, target.sample = "all")
      aipw_results[[arm_name]] <- list(
        estimate = ate["estimate"],
        std.err  = ate["std.err"],
        ci       = .make_ci(ate["estimate"], ate["std.err"])
      )

      # Per-observation influence functions for covariance computation
      scores_ipw  <- as.numeric(grf::get_scores(cf_ipw))
      scores_oadj <- as.numeric(grf::get_scores(cf_const))
      scores_aipw <- as.numeric(grf::get_scores(cf))

      # DiM influence function (block-aware propensities when blocking)
      if (!is.null(blocks_sub)) {
        p_hat_vec <- stats::ave(W_binary, blocks_sub, FUN = mean)
      } else {
        p_hat_vec <- rep(mean(W_binary), length(W_binary))
      }
      scores_dim <- W_binary / p_hat_vec * (Y_sub - mean(Y1)) -
                    (1 - W_binary) / (1 - p_hat_vec) * (Y_sub - mean(Y0))

      # 4x4 covariance matrix of ATE estimators
      score_mat <- cbind(dim = scores_dim, ipw = scores_ipw,
                         oadj = scores_oadj, aipw = scores_aipw)
      n_sub <- length(Y_sub)
      if (!is.null(clusters_sub)) {
        cluster_score_sums <- rowsum(score_mat, clusters_sub)
        n_clusters <- nrow(cluster_score_sums)
        ate_cov_list[[arm_name]] <- stats::cov(cluster_score_sums) / (n_clusters^2)
      } else {
        ate_cov_list[[arm_name]] <- stats::cov(score_mat) / n_sub
      }

      # Overlap-weighted estimates when extreme propensity scores detected
      if (has_overlap_issue) {
        ate_ipw_ow  <- grf::average_treatment_effect(cf_ipw,   target.sample = "overlap")
        ate_oadj_ow <- grf::average_treatment_effect(cf_const, target.sample = "overlap")
        ate_aipw_ow <- grf::average_treatment_effect(cf,       target.sample = "overlap")
        overlap_results[[arm_name]] <- list(
          ipw = list(
            estimate = ate_ipw_ow["estimate"],
            std.err  = ate_ipw_ow["std.err"],
            ci       = .make_ci(ate_ipw_ow["estimate"], ate_ipw_ow["std.err"])
          ),
          oadj = list(
            estimate = ate_oadj_ow["estimate"],
            std.err  = ate_oadj_ow["std.err"],
            ci       = .make_ci(ate_oadj_ow["estimate"], ate_oadj_ow["std.err"])
          ),
          aipw = list(
            estimate = ate_aipw_ow["estimate"],
            std.err  = ate_aipw_ow["std.err"],
            ci       = .make_ci(ate_aipw_ow["estimate"], ate_aipw_ow["std.err"])
          )
        )
      }
    }
  }

  # Build result object
  # For binary treatment, unwrap single-element lists for backward compatibility
  if (!multiarm) {
    result <- list(
      balance_test = balance_tests[[1]],
      dim = if (length(dim_results) > 0) dim_results[[1]] else NULL,
      ipw = if (length(ipw_results) > 0) ipw_results[[1]] else NULL,
      aipw = if (length(aipw_results) > 0) aipw_results[[1]] else NULL,
      oadj = if (length(oadj_results) > 0) oadj_results[[1]] else NULL,
      passed = unname(passed_vec[1]),
      alpha = alpha,
      cf = if (length(cf_list) > 0) cf_list[[1]] else NULL,
      ate_cov = if (length(ate_cov_list) > 0) ate_cov_list[[1]] else NULL,
      overlap_flag = overlap_flag_list[[1]],
      overlap = if (length(overlap_results) > 0) overlap_results[[1]] else NULL,
      n_extreme = n_extreme_list[[1]],
      pscores_real = pscores_real_list[[1]],
      pscores_null = pscores_null_list[[1]],
      imp.predictors = vip_list[[1]],
      n = length(W),
      n_treated = n_per_arm[[1]]$treated,
      n_control = n_per_arm[[1]]$control,
      control = control,
      arms = treatment_arms,
      multiarm = FALSE,
      clusters = clusters,
      blocks = blocks
    )
  } else {
    result <- list(
      balance_test = balance_test_joint,
      dim = if (length(dim_results) > 0) dim_results else NULL,
      ipw = if (length(ipw_results) > 0) ipw_results else NULL,
      aipw = if (length(aipw_results) > 0) aipw_results else NULL,
      oadj = if (length(oadj_results) > 0) oadj_results else NULL,
      passed = passed_joint,
      alpha = alpha,
      cf = if (length(cf_list) > 0) cf_list else NULL,
      ate_cov = if (length(ate_cov_list) > 0) ate_cov_list else NULL,
      overlap_flag = overlap_flag_list,
      overlap = if (length(overlap_results) > 0) overlap_results else NULL,
      n_extreme = n_extreme_list,
      pscores_real = pscores_real_list,
      pscores_null = pscores_null_list,
      imp.predictors = vip_list,
      n = length(W),
      n_per_arm = n_per_arm,
      control = control,
      arms = treatment_arms,
      multiarm = TRUE,
      clusters = clusters,
      blocks = blocks
    )
  }

  class(result) <- "balance"
  return(result)
}

#' @rdname balance
#' @export
print.balance <- function(x, ...) {
  LINE <- "------------------------------------------------------------\n"

  cat("\nBalance Assessment\n")
  cat(LINE)

  if (!x$multiarm) {
    status <- ifelse(x$passed, "PASS", "FAIL")
    cat(sprintf("  Control:  '%s'\n", x$control))
    if (!is.null(x$clusters) || !is.null(x$blocks)) {
      design_parts <- character(0)
      if (!is.null(x$clusters)) design_parts <- c(design_parts, sprintf("%d clusters", length(unique(x$clusters))))
      if (!is.null(x$blocks)) design_parts <- c(design_parts, sprintf("%d blocks", length(unique(x$blocks))))
      cat(sprintf("  Design:   %s\n", paste(design_parts, collapse = ", ")))
    }
    cat(sprintf("  Balance:  p = %.4f  [%s]\n", x$balance_test$pval, status))
    if (!is.null(x$dim)) {
      cat("\nTreatment Effect Estimates\n")
      cat(LINE)
      cat(sprintf("  %-28s %8.4f  (SE: %7.4f)\n", "DiM:", x$dim$estimate, x$dim$std.err))
      cat(sprintf("  %-28s %8.4f  (SE: %7.4f)\n", "IPW:", x$ipw$estimate, x$ipw$std.err))
      cat(sprintf("  %-28s %8.4f  (SE: %7.4f)\n", "Outcome-adjusted:", x$oadj$estimate, x$oadj$std.err))
      cat(sprintf("  %-28s %8.4f  (SE: %7.4f)\n", "AIPW:", x$aipw$estimate, x$aipw$std.err))
      if (isTRUE(x$overlap_flag)) {
        cat(sprintf("\n  OVERLAP WARNING: %d observations have extreme propensity scores.\n", x$n_extreme))
      }
    } else {
      cat("  No outcome Y provided: treatment effect estimates skipped.\n")
    }
  } else {
    cat(sprintf("  Control: '%s'  |  %d treatment arms\n", x$control, length(x$arms)))
    status <- ifelse(x$passed, "PASS", "FAIL")
    cat(sprintf("  Balance (joint):  p = %.4f  [%s]\n", x$balance_test$pval, status))
    cat("\n")
    has_effects <- !is.null(x$dim)
    if (has_effects) {
      for (arm in x$arms) {
        cat(sprintf("  [%s vs %s]\n", arm, x$control))
        cat(sprintf("    DiM: %8.4f   IPW: %8.4f   AIPW: %8.4f\n",
                    x$dim[[arm]]$estimate, x$ipw[[arm]]$estimate, x$aipw[[arm]]$estimate))
        if (isTRUE(x$overlap_flag[[arm]])) {
          cat(sprintf("    OVERLAP WARNING: %d extreme propensity scores\n", x$n_extreme[[arm]]))
        }
      }
    } else {
      cat("  No outcome Y provided: treatment effect estimates skipped.\n")
    }
  }

  cat("\nUse summary() for full details, plot() to visualize.\n\n")
  invisible(x)
}

#' @param object A balance result object (for summary method).
#' @rdname balance
#' @export
summary.balance <- function(object, ...) {

  SEP  <- "========================================================================\n"
  LINE <- "------------------------------------------------------------------------\n"
  has_effects <- !is.null(object$dim)

  # ── Helper: balance test block (test stats + PS diagnostics) ─────────────
  print_balance_block <- function(bt, ps_real, ps_null, passed, alpha) {
    status <- ifelse(passed, "PASS", "FAIL")
    cat(sprintf("   Classifier:          %s\n", paste(bt$class.methods, collapse = ", ")))
    cat(sprintf("   Permutations:        %d\n", bt$perm.N))
    cat(sprintf("   Test statistic:      %.4f\n", bt$teststat))
    cat(sprintf("   Null mean (SD):      %.4f (%.4f)\n",
                mean(bt$nulldist), stats::sd(bt$nulldist)))
    cat(sprintf("   P-value:             %.4f\n", bt$pval))
    cat(sprintf("   Alpha:               %.2f\n", alpha))
    cat(sprintf("   Result:              %s\n", status))
    cat("\n")
    cat("   Propensity scores (boosted regression forest):\n")
    cat(sprintf("   %-16s  %10s  %10s\n", "", "Real", "Null"))
    cat(sprintf("   %s\n", strrep("-", 40)))
    cat(sprintf("   %-16s  %10.4f  %10.4f\n", "Mean:",
                mean(ps_real), mean(ps_null)))
    cat(sprintf("   %-16s  %10.4f  %10.4f\n", "SD:",
                stats::sd(ps_real), stats::sd(ps_null)))
    cat(sprintf("   %-16s  %10.4f  %10.4f\n", "Min:",
                min(ps_real), min(ps_null)))
    cat(sprintf("   %-16s  %10.4f  %10.4f\n", "Max:",
                max(ps_real), max(ps_null)))
    cat(sprintf("   %s\n", strrep("-", 40)))
    cat(sprintf("   Diff. in means:      %.4f\n",
                mean(ps_real) - mean(ps_null)))
    cat(sprintf("   Ratio of SDs:        %.4f\n",
                stats::sd(ps_real) / stats::sd(ps_null)))
    cat("\n")
  }

  # ── Helper: estimates table rows ──────────────────────────────────────────
  # Format: label (26), estimate (9), SE (8), CI (20) = 66 chars after "   "
  EST_ROW <- "   %-26s  %9.4f  %8.4f  [%8.4f, %8.4f]\n"

  print_estimates_rows <- function(dim_res, ipw_res, aipw_res, oadj_res) {
    cat(sprintf(EST_ROW, "DiM",
                dim_res$estimate,        dim_res$std.err,
                dim_res$ci[1],           dim_res$ci[2]))
    cat(sprintf(EST_ROW, "IPW",
                ipw_res$estimate,        ipw_res$std.err,
                ipw_res$ci[1],           ipw_res$ci[2]))
    cat(sprintf(EST_ROW, "Outcome-adjusted",
                oadj_res$estimate, oadj_res$std.err,
                oadj_res$ci[1],    oadj_res$ci[2]))
    cat(sprintf(EST_ROW, "AIPW",
                aipw_res$estimate,       aipw_res$std.err,
                aipw_res$ci[1],          aipw_res$ci[2]))
  }

  # ── Helper: divergence table + dynamic interpretations ────────────────────
  # Label col = 28 to fit "AIPW (prop.) vs AIPW (dr)" etc.
  # Total: 3 + 28 + 2 + 10 + 2 + 9 + 2 + 7 + 2 + 8 = 73 chars (+ 2 sig marker)
  DIV_ROW <- "   %-28s  %10.4f  %9.4f  %7.3f  %8.4f%s\n"

  # Estimator index map for covariance matrix: dim=1, ipw=2, oadj=3, aipw=4
  print_divergence_rows <- function(dim_res, ipw_res, aipw_res, oadj_res, alpha, Sigma) {
    pairs <- list(
      list(key = "dim_ipw", label = "DiM vs IPW",
           e1 = dim_res$estimate,        e2 = ipw_res$estimate,
           idx1 = 1L, idx2 = 2L),
      list(key = "dim_oa",  label = "DiM vs Outcome-adj.",
           e1 = dim_res$estimate,        e2 = oadj_res$estimate,
           idx1 = 1L, idx2 = 3L),
      list(key = "dim_dr",  label = "DiM vs AIPW",
           e1 = dim_res$estimate,        e2 = aipw_res$estimate,
           idx1 = 1L, idx2 = 4L),
      list(key = "ipw_dr",  label = "IPW vs AIPW",
           e1 = ipw_res$estimate,        e2 = aipw_res$estimate,
           idx1 = 2L, idx2 = 4L)
    )

    # Compute stats for all pairs using the covariance matrix
    res <- lapply(pairs, function(p) {
      d  <- p$e1 - p$e2
      se <- sqrt(Sigma[p$idx1, p$idx1] + Sigma[p$idx2, p$idx2] - 2 * Sigma[p$idx1, p$idx2])
      z  <- d / se
      pv <- 2 * stats::pnorm(-abs(z))
      list(key = p$key, label = p$label,
           diff = d, se = se, z = z, pval = pv, sig = pv < alpha)
    })
    r <- stats::setNames(res, vapply(res, `[[`, character(1L), "key"))

    # Print table
    any_sig <- any(vapply(res, `[[`, logical(1L), "sig"))
    for (rv in res) {
      cat(sprintf(DIV_ROW, rv$label, rv$diff, rv$se, rv$z, rv$pval,
                  if (rv$sig) " *" else "  "))
    }
    if (any_sig) cat(sprintf("   * p < %.2f\n", alpha))
    cat("\n")

    # Robustness summary
    if (!any_sig) {
      cat("   All four estimators agree closely, indicating the ATE estimate is\n")
      cat("   robust to the choice of nuisance model specification.\n\n")
      return(invisible(NULL))
    }

    cat("   Significant divergences:\n\n")

    # a) DiM vs IPW
    if (r$dim_ipw$sig) {
      cat(sprintf(paste0(
        "   DiM vs IPW: The IPW estimate differs from the unadjusted DiM\n",
        "   by %.4f units (z = %.3f, p = %.4f), indicating that propensity\n",
        "   reweighting accounts for this difference.\n\n"
      ), abs(r$dim_ipw$diff), r$dim_ipw$z, r$dim_ipw$pval))
    }

    # b) DiM vs Outcome-adjusted
    if (r$dim_oa$sig) {
      cat(sprintf(paste0(
        "   DiM vs Outcome-adjusted: The outcome-adjusted estimate differs from\n",
        "   the unadjusted DiM by %.4f units (z = %.3f, p = %.4f), indicating\n",
        "   that outcome regression adjustment accounts for this difference.\n\n"
      ), abs(r$dim_oa$diff), r$dim_oa$z, r$dim_oa$pval))
    }

    # c) DiM vs AIPW
    if (r$dim_dr$sig) {
      cat(sprintf(paste0(
        "   DiM vs AIPW: The AIPW estimate differs from the unadjusted DiM\n",
        "   by %.4f units (z = %.3f, p = %.4f), reflecting the combined\n",
        "   effect of propensity and outcome adjustment.\n\n"
      ), abs(r$dim_dr$diff), r$dim_dr$z, r$dim_dr$pval))
    }

    # d) IPW vs AIPW
    if (r$ipw_dr$sig) {
      cat(sprintf(paste0(
        "   IPW vs AIPW: Adding outcome regression to propensity reweighting\n",
        "   changes the estimate by %.4f units (z = %.3f, p = %.4f),\n",
        "   indicating that outcome modeling captures additional\n",
        "   covariate-outcome associations beyond reweighting.\n\n"
      ), abs(r$ipw_dr$diff), r$ipw_dr$z, r$ipw_dr$pval))
    }
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PART I: COVARIATE BALANCE ASSESSMENT
  # ══════════════════════════════════════════════════════════════════════════
  cat("\n")
  cat(SEP)
  cat("                   COVARIATE BALANCE ASSESSMENT                        \n")
  cat(SEP)
  cat("\n")

  # Section 1: Sample
  cat("1. SAMPLE\n")
  cat(LINE)
  cat(sprintf("   Observations:    %d\n", object$n))
  if (!is.null(object$clusters))
    cat(sprintf("   Clusters:        %d\n", length(unique(object$clusters))))
  if (!is.null(object$blocks))
    cat(sprintf("   Blocks:          %d\n", length(unique(object$blocks))))
  if (!object$multiarm) {
    cat(sprintf("   Control ('%s'):  %d (%.1f%%)\n",
                object$control, object$n_control,
                100 * object$n_control / object$n))
    cat(sprintf("   Treatment:       %d (%.1f%%)\n",
                object$n_treated, 100 * object$n_treated / object$n))
  } else {
    n_ctrl <- object$n_per_arm[[1]]$control
    cat(sprintf("   Control ('%s'):  %d (%.1f%%)\n",
                object$control, n_ctrl, 100 * n_ctrl / object$n))
    for (arm in object$arms) {
      n_arm <- object$n_per_arm[[arm]]$treated
      cat(sprintf("   Arm ('%s'):  %d (%.1f%%)\n",
                  arm, n_arm, 100 * n_arm / object$n))
    }
  }
  cat("\n")

  if (!object$multiarm) {
    # ── Binary: section 2 = balance test, section 3 = interpretation ─────
    cat("2. CLASSIFICATION PERMUTATION TEST\n")
    cat(LINE)
    print_balance_block(object$balance_test, object$pscores_real, object$pscores_null,
                        object$passed, object$alpha)

    cat("3. INTERPRETATION\n")
    cat(LINE)
    if (object$passed) {
      cat("   The classification permutation test does not reject the null\n")
      cat("   hypothesis that treatment and control groups are drawn from the\n")
      cat("   same covariate distribution. The classifier cannot distinguish\n")
      cat("   between groups better than chance.\n")
    } else {
      cat("   The classification permutation test rejects the null hypothesis,\n")
      cat("   indicating that treatment and control groups differ in their\n")
      cat("   covariate distributions. The classifier can distinguish between\n")
      cat("   groups better than random chance.\n")
    }
    cat("\n")
    next_sec <- 4L

  } else {
    # ── Multi-arm: section 2 = joint K-class test, section 3 = per-arm PS diagnostics ──
    bt <- object$balance_test
    cat("2. CLASSIFICATION PERMUTATION TEST (JOINT)\n")
    cat(LINE)
    status <- ifelse(object$passed, "PASS", "FAIL")
    cat(sprintf("   Type:                Joint %d-class test\n", length(object$arms) + 1L))
    cat(sprintf("   Classifier:          %s\n", paste(bt$class.methods, collapse = ", ")))
    cat(sprintf("   Permutations:        %d\n", bt$perm.N))
    cat(sprintf("   Test statistic:      %.4f\n", bt$teststat))
    cat(sprintf("   Null mean (SD):      %.4f (%.4f)\n",
                mean(bt$nulldist), stats::sd(bt$nulldist)))
    cat(sprintf("   P-value:             %.4f\n", bt$pval))
    cat(sprintf("   Alpha:               %.2f\n", object$alpha))
    cat(sprintf("   Result:              %s\n", status))
    cat("\n")

    sec <- 3L
    for (arm in object$arms) {
      cat(sprintf("%d. PROPENSITY DIAGNOSTICS: %s vs %s\n", sec, arm, object$control))
      cat(LINE)
      ps_real <- object$pscores_real[[arm]]
      ps_null <- object$pscores_null[[arm]]
      cat("   Propensity scores (boosted regression forest):\n")
      cat(sprintf("   %-16s  %10s  %10s\n", "", "Real", "Null"))
      cat(sprintf("   %s\n", strrep("-", 40)))
      cat(sprintf("   %-16s  %10.4f  %10.4f\n", "Mean:", mean(ps_real), mean(ps_null)))
      cat(sprintf("   %-16s  %10.4f  %10.4f\n", "SD:", stats::sd(ps_real), stats::sd(ps_null)))
      cat(sprintf("   %-16s  %10.4f  %10.4f\n", "Min:", min(ps_real), min(ps_null)))
      cat(sprintf("   %-16s  %10.4f  %10.4f\n", "Max:", max(ps_real), max(ps_null)))
      cat(sprintf("   %s\n", strrep("-", 40)))
      cat(sprintf("   Diff. in means:      %.4f\n", mean(ps_real) - mean(ps_null)))
      cat(sprintf("   Ratio of SDs:        %.4f\n", stats::sd(ps_real) / stats::sd(ps_null)))
      cat("\n")
      sec <- sec + 1L
    }

    cat(sprintf("%d. INTERPRETATION\n", sec))
    cat(LINE)
    if (object$passed) {
      cat("   The joint classification permutation test does not reject the null\n")
      cat("   hypothesis that all treatment groups are drawn from the same\n")
      cat("   covariate distribution. The classifier cannot distinguish between\n")
      cat("   groups better than chance.\n")
    } else {
      cat("   The joint classification permutation test rejects the null\n")
      cat("   hypothesis, indicating that at least one treatment group differs\n")
      cat("   in its covariate distribution. The classifier can distinguish\n")
      cat("   between groups better than random chance.\n")
    }
    cat("\n")
    next_sec <- sec + 1L
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PART II: TREATMENT EFFECT ESTIMATION
  # ══════════════════════════════════════════════════════════════════════════
  if (!has_effects) {
    return(invisible(object))
  }

  cat(SEP)
  cat("                   TREATMENT EFFECT ESTIMATION                         \n")
  cat(SEP)
  cat("\n")

  # Shared table header/separator strings
  # EST_ROW: 3 + 26 + 2 + 9 + 2 + 8 + 2 + [1+8+2+8+1] = 72 chars
  EST_HEAD <- sprintf("   %-26s  %9s  %8s  %20s\n",
                      "Estimator", "Estimate", "SE", "95% CI")
  EST_SEP  <- sprintf("   %s\n", strrep("-", 69))
  # DIV_ROW: 3 + 28 + 2 + 10 + 2 + 9 + 2 + 7 + 2 + 8 = 73 chars (+ 2 for sig marker)
  DIV_HEAD <- sprintf("   %-28s  %10s  %9s  %7s  %8s\n",
                      "Comparison", "Difference", "SE(diff)", "z-stat", "p-value")
  DIV_SEP  <- sprintf("   %s\n", strrep("-", 70))

  # ── Helper: overlap warning block ────────────────────────────────────────
  OW_ROW <- "   %-26s  %9.4f  %8.4f  [%8.4f, %8.4f]\n"
  print_overlap_block <- function(n_ext, ov) {
    cat(sprintf("\n   OVERLAP WARNING: %d observations have extreme propensity scores\n", n_ext))
    cat("   (< 0.05 or > 0.95). Overlap-weighted estimates down-weight these:\n\n")
    cat(EST_HEAD)
    cat(EST_SEP)
    cat(sprintf(OW_ROW, "IPW (OW)",
                ov$ipw$estimate, ov$ipw$std.err, ov$ipw$ci[1], ov$ipw$ci[2]))
    cat(sprintf(OW_ROW, "Outcome-adj. (OW)",
                ov$oadj$estimate, ov$oadj$std.err, ov$oadj$ci[1], ov$oadj$ci[2]))
    cat(sprintf(OW_ROW, "AIPW (OW)",
                ov$aipw$estimate, ov$aipw$std.err, ov$aipw$ci[1], ov$aipw$ci[2]))
    cat("\n")
  }

  if (!object$multiarm) {
    # Section N: Estimates
    cat(sprintf("%d. ESTIMATES\n", next_sec))
    cat(LINE)
    cat(EST_HEAD)
    cat(EST_SEP)
    print_estimates_rows(object$dim, object$ipw, object$aipw, object$oadj)
    if (isTRUE(object$overlap_flag) && !is.null(object$overlap)) {
      print_overlap_block(object$n_extreme, object$overlap)
    }
    cat("\n")

    # Section N+1: Divergence tests
    cat(sprintf("%d. ESTIMATOR DIVERGENCE TESTS\n", next_sec + 1L))
    cat(LINE)
    cat(DIV_HEAD)
    cat(DIV_SEP)
    print_divergence_rows(object$dim, object$ipw, object$aipw, object$oadj, object$alpha, object$ate_cov)

    next_sec <- next_sec + 2L

  } else {
    # Estimates: one section, arm sub-headers inside the table
    cat(sprintf("%d. ESTIMATES\n", next_sec))
    cat(LINE)
    cat(EST_HEAD)
    cat(EST_SEP)
    for (arm in object$arms) {
      cat(sprintf("   --- %s vs %s ---\n", arm, object$control))
      print_estimates_rows(object$dim[[arm]], object$ipw[[arm]],
                           object$aipw[[arm]], object$oadj[[arm]])
      if (isTRUE(object$overlap_flag[[arm]]) && !is.null(object$overlap[[arm]])) {
        print_overlap_block(object$n_extreme[[arm]], object$overlap[[arm]])
      }
    }
    cat("\n")

    # Divergence tests: one section, arm sub-headers inside
    cat(sprintf("%d. ESTIMATOR DIVERGENCE TESTS\n", next_sec + 1L))
    cat(LINE)
    for (arm in object$arms) {
      cat(sprintf("   --- %s vs %s ---\n", arm, object$control))
      cat(DIV_HEAD)
      cat(DIV_SEP)
      print_divergence_rows(object$dim[[arm]], object$ipw[[arm]],
                            object$aipw[[arm]], object$oadj[[arm]], object$alpha,
                            object$ate_cov[[arm]])
    }

    next_sec <- next_sec + 2L
  }

  # Final section: Guide
  cat(sprintf("%d. ESTIMATOR GUIDE\n", next_sec))
  cat(LINE)
  cat("   All adjusted estimators use grf::causal_forest's AIPW framework.\n")
  cat("   They differ in which nuisance models (propensity and/or outcome)\n")
  cat("   are estimated vs. held at uninformative constants.\n")
  cat("\n")
  cat("   We present four ATE estimates as a robustness decomposition.\n")
  cat("   Agreement across estimators signals robustness to modeling choices.\n")
  cat("   Divergence reveals which adjustment component (propensity vs.\n")
  cat("   outcome) most affects the estimate and warrants investigation.\n")
  cat("\n")
  cat("   Nuisance models are fit using honest, tuned boosted regression\n")
  cat("   forests (grf::boosted_regression_forest; honesty = TRUE,\n")
  cat("   tune.parameters = \"all\"); all predictions are out-of-bag.\n")
  cat("   Treatment effects are estimated via grf::causal_forest\n")
  cat("   and grf::average_treatment_effect (target = \"all\").\n")
  cat("   Standard errors use the infinitesimal jackknife (IJ) influence\n")
  cat("   function. DiM SEs use Neyman's separate-variance (Welch) formula.\n")
  if (!is.null(object$clusters)) {
    cat("   Cluster-robust variance estimators are used throughout (DiM SEs\n")
    cat("   use cluster means; causal forest SEs account for clustering;\n")
    cat("   divergence test covariance aggregates scores to cluster level).\n")
  }
  cat("   All confidence intervals use a normal approximation.\n")
  cat("\n")
  cat("   DiM (difference-in-means)\n")
  cat("     E[Y|W=1] - E[Y|W=0]. No covariate adjustment.\n")
  cat("\n")
  cat("   IPW (inverse propensity weighted)\n")
  cat("     W.hat = boosted-RF propensity; Y.hat = mean(Y) (constant).\n")
  cat("     Isolates the effect of propensity score reweighting.\n")
  cat("\n")
  cat("   Outcome-adjusted (regression adjustment)\n")
  cat("     W.hat = mean(W) (constant); Y.hat = boosted-RF outcome predictions.\n")
  cat("     Isolates the effect of outcome regression adjustment.\n")
  cat("\n")
  cat("   AIPW (augmented IPW / doubly robust)\n")
  cat("     Both nuisance models estimated. Consistent if either is correct.\n")
  cat("\n")

  # Check if any overlap warnings were triggered
  has_overlap <- if (!object$multiarm) {
    isTRUE(object$overlap_flag) && !is.null(object$overlap)
  } else {
    !is.null(object$overlap) &&
      any(vapply(object$arms, function(a) isTRUE(object$overlap_flag[[a]]), logical(1L)))
  }
  if (has_overlap) {
    cat("   Overlap-Weighted (OW) Estimates\n")
    cat("     When propensity scores approach 0 or 1, AIPW estimators become\n")
    cat("     unstable. Overlap-weighted estimates (Li, Morgan & Zaslavsky, 2018)\n")
    cat("     target the ATE for the overlap population by down-weighting units\n")
    cat("     with extreme propensity scores. Computed via\n")
    cat("     grf::average_treatment_effect(target.sample = \"overlap\").\n")
    cat("\n")
  }

  invisible(object)
}

#' Plot method for balance objects
#'
#' @param x A balance result object.
#' @param which Character vector specifying which plots to create. Options are "pscores", "null_dist", "effects", or "all".
#' @param combined Logical. If TRUE, displays all three plots in a combined panel. Default is TRUE.
#' @param breaks Number of breaks for histograms. Default is 25.
#' @param ... Additional arguments (currently unused).
#' @rdname balance
#' @export
plot.balance <- function(x, which = "all", combined = TRUE, breaks = 25, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for plotting.", call. = FALSE)
  }

  # Color palette (consistent across plots)
  col_null <- "#4575B4"      # Blue for null/reference
  col_real <- "#D73027"      # Red-orange for real/observed
  col_pass <- "#1A9850"      # Green for pass
  col_fail <- "#D73027"      # Red for fail

  # Dynamic x-axis label from metric name
  metric_labels <- c("probability" = "Probability", "rate" = "Classification Rate",
                      "mse" = "MSE", "logscore" = "Log Score")
  mn <- x$balance_test$metric_name
  stat_label <- if (!is.null(mn) && mn %in% names(metric_labels))
    paste("Test Statistic -", metric_labels[mn]) else "Test Statistic"

  # Shared theme (matching pdp style)
  g_theme <- function() {
    ggplot2::theme(
      text = ggplot2::element_text(size = 12, family = "serif"),
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
      panel.background = ggplot2::element_rect(fill = "#f0f0f0"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(color = "grey30", size = 10),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, color = "grey40"),
      legend.position = c(0.12, 0.88),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      panel.spacing = ggplot2::unit(1, "lines"),
      complete = TRUE
    )
  }

  has_effects <- !is.null(x$dim) && !is.null(x$ipw) && !is.null(x$aipw) && !is.null(x$oadj)

  if (length(which) == 1 && which == "all") {
    which <- if (has_effects) c("pscores", "null_dist", "effects") else c("pscores", "null_dist")
  }
  if ("effects" %in% which && !has_effects) {
    warning("Outcome Y not provided: skipping treatment effect plot.", call. = FALSE)
    which <- setdiff(which, "effects")
  }

  plots <- list()

  # ============================================================================
  # BINARY TREATMENT PLOTS
  # ============================================================================
  if (!x$multiarm) {

    # Panel A: Propensity score distributions
    if ("pscores" %in% which && !is.null(x$pscores_real) && !is.null(x$pscores_null)) {
      plot_df <- data.frame(
        var = factor(c(rep("Real", length(x$pscores_real)), rep("Null", length(x$pscores_null))),
                     levels = c("Null", "Real")),
        val = c(x$pscores_real, x$pscores_null)
      )

      plots$pscores <- ggplot2::ggplot(plot_df, ggplot2::aes(x = val, fill = var)) +
        ggdist::stat_histinterval(
          slab_color = "gray70",
          outline_bars = TRUE,
          alpha = 0.75,
          point_alpha = 0,
          slab_linewidth = 0.5,
          breaks = breaks,
          interval_alpha = 0
        ) +
        ggplot2::geom_vline(xintercept = mean(x$pscores_real), color = "darkorange1", linetype = "dotdash", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = mean(x$pscores_null), color = "dodgerblue1", linetype = "dotdash", linewidth = 0.5) +
        g_theme() +
        ggplot2::labs(
          title = "A. Propensity Score Distributions",
          x = "Treatment Propensity Scores",
          y = "Density",
          caption = expression(italic("Note: Dotted lines represent mean values for the null and real treatment propensity distributions."))
        ) +
        ggplot2::scale_fill_manual(values = c("dodgerblue1", "darkorange1"), name = "") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "")) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        if (max(plot_df$val) > 1 || min(plot_df$val) < 0) {
          ggplot2::scale_x_continuous(expand = c(0, 0))
        } else {
          ggplot2::scale_x_continuous(limits = c(0, 1.01), expand = c(0, 0))
        }
    }

    # Panel B: Classification permutation test null distribution
    if ("null_dist" %in% which) {
      test_pass <- x$balance_test$pval > x$alpha
      color_select <- ifelse(test_pass, col_pass, col_fail)
      status_text <- ifelse(test_pass, "Pass", "Fail")
      testval <- x$balance_test$teststat
      null_df <- data.frame(val = x$balance_test$nulldist)

      x_min <- floor(min(x$balance_test$nulldist) * 10) / 10
      x_max <- ceiling(max(c(x$balance_test$nulldist, testval)) * 10) / 10

      plots$null_dist <- ggplot2::ggplot(null_df, ggplot2::aes(x = val)) +
        ggdist::stat_histinterval(
          fill = "grey70",
          slab_color = "gray50",
          outline_bars = TRUE,
          alpha = 0.75,
          point_alpha = 0,
          slab_linewidth = 0.5,
          breaks = breaks,
          interval_alpha = 0
        ) +
        ggplot2::geom_vline(xintercept = testval, color = color_select, linewidth = 1.5) +
        ggplot2::annotate("text", x = testval, y = 0.95,
                          label = sprintf("p = %.3f (%s)", x$balance_test$pval, status_text),
                          hjust = -0.1, size = 3.5, fontface = "bold", color = color_select) +
        g_theme() +
        ggplot2::labs(
          title = "B. Classification Permutation Test",
          x = stat_label,
          y = "Density",
          caption = sprintf("Note: Classifier = %s, alpha = %.2f.",
                            paste(x$balance_test$class.methods, collapse = ", "), x$alpha)
        ) +
        ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0))
    }

    # Panel C: Treatment effect estimates with ggdist
    if ("effects" %in% which && has_effects) {
      # Colors for estimators
      col_dim        <- "#506D27"
      col_ipw        <- "#F9E79F"
      col_oadj       <- "#B2C0DA"
      col_aipw       <- "#9C4643"

      has_ow <- isTRUE(x$overlap_flag) && !is.null(x$overlap)

      # Labels matching summary output
      est_labels <- c("DiM", "IPW", "Outcome-\nadjusted", "AIPW")

      # Build data frame with weight type for shape mapping
      effect_df <- data.frame(
        estimator = factor(est_labels, levels = est_labels),
        estimate  = c(x$dim$estimate, x$ipw$estimate, x$oadj$estimate, x$aipw$estimate),
        se        = c(x$dim$std.err, x$ipw$std.err, x$oadj$std.err, x$aipw$std.err),
        weight    = "ATE"
      )

      if (has_ow) {
        ow_df <- data.frame(
          estimator = factor(est_labels[2:4], levels = est_labels),
          estimate  = c(x$overlap$ipw$estimate, x$overlap$oadj$estimate, x$overlap$aipw$estimate),
          se        = c(x$overlap$ipw$std.err, x$overlap$oadj$std.err, x$overlap$aipw$std.err),
          weight    = "ATO"
        )
        effect_df <- rbind(effect_df, ow_df)
      }
      effect_df$weight <- factor(effect_df$weight, levels = c("ATE", "ATO"))

      # Create distributional data for ggdist
      effect_df$.dist <- lapply(seq_len(nrow(effect_df)), function(i) {
        distributional::dist_normal(effect_df$estimate[i], effect_df$se[i])
      })
      effect_df$.dist <- do.call(c, effect_df$.dist)

      ow_caption <- "Note: Point estimates with 80%, 95%, 99% CIs. Intervals based on normal approximation."
      if (has_ow) ow_caption <- paste0(ow_caption, "\nOverlap-weighted estimates down-weight extreme propensity scores.")

      p <- ggplot2::ggplot(effect_df, ggplot2::aes(
        x = estimator, ydist = .dist, color = estimator, fill = estimator,
        shape = weight, group = interaction(estimator, weight)
      )) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        ggdist::stat_pointinterval(
          .width = c(0.80, 0.95, 0.99),
          point_size = 3,
          point_colour = "black",
          stroke = 0.6,
          interval_size_range = c(0.8, 2.5),
          interval_alpha = 0.8,
          position = if (has_ow) ggplot2::position_dodge(width = 0.5) else "identity"
        ) +
        ggplot2::scale_color_manual(values = c(col_dim, col_ipw, col_oadj, col_aipw), guide = "none") +
        ggplot2::scale_fill_manual(values = c(col_dim, col_ipw, col_oadj, col_aipw), guide = "none") +
        ggplot2::guides(color = "none", fill = "none") +
        g_theme() +
        ggplot2::labs(
          title = "C. Treatment Effect Estimates",
          x = NULL,
          y = "Effect Estimate",
          caption = ow_caption
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(size = 9, color = "black"),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank()
        )

      if (has_ow) {
        p <- p +
          ggplot2::scale_shape_manual(values = c("ATE" = 21, "ATO" = 24), name = "Weighting") +
          ggplot2::theme(
            legend.position = c(0.12, 0.88),
            legend.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank()
          )
      } else {
        p <- p +
          ggplot2::scale_shape_manual(values = c("ATE" = 21), guide = "none") +
          ggplot2::theme(legend.position = "none")
      }

      plots$effects <- p
    }

  } else {
    # ============================================================================
    # MULTI-ARM TREATMENT PLOTS (with facets)
    # ============================================================================

    # Panel A: Propensity score distributions (faceted by arm)
    if ("pscores" %in% which && !is.null(x$pscores_real)) {
      # Build combined data frame for all arms
      plot_df_list <- lapply(x$arms, function(arm) {
        data.frame(
          arm = sprintf("%s vs %s", arm, x$control),
          var = factor(c(rep("Real", length(x$pscores_real[[arm]])),
                         rep("Null", length(x$pscores_null[[arm]]))),
                       levels = c("Null", "Real")),
          val = c(x$pscores_real[[arm]], x$pscores_null[[arm]])
        )
      })
      plot_df <- do.call(rbind, plot_df_list)
      plot_df$arm <- factor(plot_df$arm, levels = unique(plot_df$arm))

      # Compute means per arm for vertical lines
      mean_df <- do.call(rbind, lapply(x$arms, function(arm) {
        data.frame(
          arm = sprintf("%s vs %s", arm, x$control),
          mean_real = mean(x$pscores_real[[arm]]),
          mean_null = mean(x$pscores_null[[arm]])
        )
      }))
      mean_df$arm <- factor(mean_df$arm, levels = levels(plot_df$arm))

      plots$pscores <- ggplot2::ggplot(plot_df, ggplot2::aes(x = val, fill = var)) +
        ggdist::stat_histinterval(
          slab_color = "gray70",
          outline_bars = TRUE,
          alpha = 0.75,
          point_alpha = 0,
          slab_linewidth = 0.5,
          breaks = breaks,
          interval_alpha = 0
        ) +
        ggplot2::geom_vline(data = mean_df, ggplot2::aes(xintercept = mean_real),
                            color = "darkorange1", linetype = "dotdash", linewidth = 0.5) +
        ggplot2::geom_vline(data = mean_df, ggplot2::aes(xintercept = mean_null),
                            color = "dodgerblue1", linetype = "dotdash", linewidth = 0.5) +
        ggplot2::facet_wrap(~ arm, scales = "free_x") +
        g_theme() +
        ggplot2::labs(
          title = "A. Propensity Score Distributions by Treatment Arm",
          x = "Treatment Propensity Scores",
          y = "Density",
          caption = expression(italic("Note: Dotted lines represent mean values for the null and real treatment propensity distributions."))
        ) +
        ggplot2::scale_fill_manual(values = c("dodgerblue1", "darkorange1"), name = "") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "")) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::theme(legend.position = "bottom")
    }

    # Panel B: Classification permutation test null distribution (single joint test)
    if ("null_dist" %in% which) {
      test_pass <- x$balance_test$pval > x$alpha
      color_select <- ifelse(test_pass, col_pass, col_fail)
      status_text <- ifelse(test_pass, "Pass", "Fail")
      testval <- x$balance_test$teststat
      null_df <- data.frame(val = x$balance_test$nulldist)

      x_min <- floor(min(x$balance_test$nulldist) * 10) / 10
      x_max <- ceiling(max(c(x$balance_test$nulldist, testval)) * 10) / 10

      plots$null_dist <- ggplot2::ggplot(null_df, ggplot2::aes(x = val)) +
        ggdist::stat_histinterval(
          fill = "grey70",
          slab_color = "gray50",
          outline_bars = TRUE,
          alpha = 0.75,
          point_alpha = 0,
          slab_linewidth = 0.5,
          breaks = breaks,
          interval_alpha = 0
        ) +
        ggplot2::geom_vline(xintercept = testval, color = color_select, linewidth = 1.5) +
        ggplot2::annotate("text", x = testval, y = 0.95,
                          label = sprintf("p = %.3f (%s)", x$balance_test$pval, status_text),
                          hjust = -0.1, size = 3.5, fontface = "bold", color = color_select) +
        g_theme() +
        ggplot2::labs(
          title = "B. Joint Classification Permutation Test",
          x = stat_label,
          y = "Density",
          caption = sprintf("Note: Joint %d-class test. Classifier = %s, alpha = %.2f.",
                            length(x$arms) + 1L,
                            paste(x$balance_test$class.methods, collapse = ", "), x$alpha)
        ) +
        ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0))
    }

    # Panel C: Treatment effect estimates (faceted by arm)
    if ("effects" %in% which && has_effects) {
      # Colors for estimators
      col_dim        <- "#506D27"
      col_ipw        <- "#F9E79F"
      col_oadj       <- "#B2C0DA"
      col_aipw       <- "#9C4643"

      any_overlap <- any(vapply(x$arms, function(arm) isTRUE(x$overlap_flag[[arm]]), logical(1L)))

      # Labels matching summary output
      est_labels <- c("DiM", "IPW", "Outcome-\nadjusted", "AIPW")

      # Build combined data frame for all arms with weight type
      effect_df_list <- lapply(x$arms, function(arm) {
        base_df <- data.frame(
          arm       = sprintf("%s vs %s", arm, x$control),
          estimator = factor(est_labels, levels = est_labels),
          estimate  = c(x$dim[[arm]]$estimate, x$ipw[[arm]]$estimate, x$oadj[[arm]]$estimate, x$aipw[[arm]]$estimate),
          se        = c(x$dim[[arm]]$std.err, x$ipw[[arm]]$std.err, x$oadj[[arm]]$std.err, x$aipw[[arm]]$std.err),
          weight    = "ATE"
        )
        if (isTRUE(x$overlap_flag[[arm]]) && !is.null(x$overlap) && !is.null(x$overlap[[arm]])) {
          ow <- x$overlap[[arm]]
          ow_df <- data.frame(
            arm       = sprintf("%s vs %s", arm, x$control),
            estimator = factor(est_labels[2:4], levels = est_labels),
            estimate  = c(ow$ipw$estimate, ow$oadj$estimate, ow$aipw$estimate),
            se        = c(ow$ipw$std.err, ow$oadj$std.err, ow$aipw$std.err),
            weight    = "ATO"
          )
          base_df <- rbind(base_df, ow_df)
        }
        base_df
      })
      effect_df <- do.call(rbind, effect_df_list)
      effect_df$arm <- factor(effect_df$arm, levels = unique(effect_df$arm))
      effect_df$weight <- factor(effect_df$weight, levels = c("ATE", "ATO"))

      # Create distributional data for ggdist
      effect_df$.dist <- lapply(seq_len(nrow(effect_df)), function(i) {
        distributional::dist_normal(effect_df$estimate[i], effect_df$se[i])
      })
      effect_df$.dist <- do.call(c, effect_df$.dist)

      ow_caption <- "Note: Point estimates with 80%, 95%, 99% CIs. Intervals based on normal approximation."
      if (any_overlap) ow_caption <- paste0(ow_caption, "\nOverlap-weighted estimates down-weight extreme propensity scores.")

      p <- ggplot2::ggplot(effect_df, ggplot2::aes(
        x = estimator, ydist = .dist, color = estimator, fill = estimator,
        shape = weight, group = interaction(estimator, weight)
      )) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        ggdist::stat_pointinterval(
          .width = c(0.80, 0.95, 0.99),
          point_size = 3,
          point_colour = "black",
          stroke = 0.6,
          interval_size_range = c(0.8, 2.5),
          interval_alpha = 0.8,
          position = if (any_overlap) ggplot2::position_dodge(width = 0.5) else "identity"
        ) +
        ggplot2::scale_color_manual(values = c(col_dim, col_ipw, col_oadj, col_aipw), guide = "none") +
        ggplot2::scale_fill_manual(values = c(col_dim, col_ipw, col_oadj, col_aipw), guide = "none") +
        ggplot2::guides(color = "none", fill = "none") +
        ggplot2::facet_wrap(~ arm) +
        g_theme() +
        ggplot2::labs(
          title = "C. Treatment Effect Estimates by Arm",
          x = NULL,
          y = "Effect Estimate",
          caption = ow_caption
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(size = 8, color = "black", angle = 45, hjust = 1),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank()
        )

      if (any_overlap) {
        p <- p +
          ggplot2::scale_shape_manual(values = c("ATE" = 21, "ATO" = 24), name = "Weighting") +
          ggplot2::theme(
            legend.position = "bottom",
            legend.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank()
          )
      } else {
        p <- p +
          ggplot2::scale_shape_manual(values = c("ATE" = 21), guide = "none") +
          ggplot2::theme(legend.position = "none")
      }

      plots$effects <- p
    }
  }

  # Combined multi-panel display
  if (combined && length(plots) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required for combined plots.", call. = FALSE)
    }
    # For multi-arm, stack vertically since each plot is already wide with facets
    if (x$multiarm) {
      return(patchwork::wrap_plots(plotlist = unname(plots), ncol = 1))
    } else {
      return(patchwork::wrap_plots(plotlist = unname(plots)))
    }
  }

  if (length(plots) == 1) {
    return(plots[[1]])
  }

  return(plots)
}

