# Suppress R CMD check notes for non-standard evaluation
utils::globalVariables(c("estimator", "estimate", "var", "val", ".dist", "arm", "comparison",
                         "mean_real", "mean_null", "testval", "pval", "status_text", "color"))

#' Complete Balance Assessment and Treatment Effect Estimation
#'
#' @description
#' A unified function for covariate balance assessment and treatment effect estimation.
#' Combines a formal balance test (via classification permutation test), visual diagnostics
#' (propensity score distributions), and treatment effect estimates using both difference-in-means
#' and augmented inverse propensity weighting (AIPW).
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
#' @param class.method Classification method for balance test. Default is "ferns".
#' @param seed Random seed for reproducibility. Default is 1995.
#' @param control Optional. The value in \code{W} to use as the control group. If \code{NULL} (default),
#'   the first factor level is used as control. A message is displayed indicating the control assumption.
#' @param fastcpt.args A named list of additional arguments to pass to \code{\link{fastcpt}}. For example, \code{fastcpt.args = list(parallel = TRUE, leaveout = 0.2)}.
#'
#' @return A list of class "balance" containing:
#' \item{balance_test}{Results from fastcpt including p-value and propensity scores. For multi-arm, a named list with one entry per treatment arm.}
#' \item{dim}{Difference-in-means estimate with standard error and confidence interval (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{aipw}{Doubly robust estimate from causal forest (with propensity weighting) with standard error and CI (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{aipw_const}{Outcome-adjusted estimate from causal forest (no propensity weighting) with SE and CI (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{passed}{Logical indicating whether the balance test passed. For multi-arm, a named logical vector.}
#' \item{alpha}{The significance level used.}
#' \item{cf}{The fitted causal_forest object(s) for advanced users (only if \code{Y} is provided). For multi-arm, a named list.}
#' \item{control}{The control level used.}
#' \item{arms}{Character vector of treatment arm names (excluding control).}
#' \item{multiarm}{Logical indicating whether this is a multi-arm analysis.}
#'
#' @examples
#' \dontrun{
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
                    control = NULL, fastcpt.args = list()) {
  
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
  if (!is.null(Y)) {
    if (!is.numeric(Y)) stop("Y must be a numeric vector.", call. = FALSE)
    if (length(Y) != length(W)) {
      stop("Y and W must have the same number of observations.", call. = FALSE)
    }
  }
  
  # Determine control level
  W_factor <- as.factor(W)
  all_levels <- levels(W_factor)
  
  if (is.null(control)) {
    control <- all_levels[1]
    if (n_arms > 2) {
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
    message(sprintf("Using '%s' as control group.", control))
  }
  
  # Treatment arms (excluding control)
  treatment_arms <- setdiff(all_levels, control)
  multiarm <- n_arms > 2
  
  # Convert X to matrix, handling factors for grf
  X_df <- as.data.frame(X)
  has_factors <- any(sapply(X_df, is.factor) | sapply(X_df, is.character))
  
  if (has_factors) {
    # Expand factors to dummies for grf
    X_matrix <- stats::model.matrix(~ . - 1, data = X_df)
    # Check for high-cardinality factors
    n_cols <- ncol(X_matrix)
    if (n_cols > ncol(X_df) * 5) {
      warning("Factor expansion created many columns. Consider collapsing high-cardinality factors.", 
              call. = FALSE)
    }
  } else {
    X_matrix <- as.matrix(X_df)
  }
  
  # Set seed
  set.seed(seed)
  
  # Initialize storage for results
  balance_tests <- list()
  pscores_real_list <- list()
  pscores_null_list <- list()
  dim_results <- list()
  aipw_results <- list()
  aipw_const_results <- list()
  cf_list <- list()
  passed_vec <- c()
  n_per_arm <- list()
  
  # Loop over each treatment arm (pairwise vs control)
  for (arm in treatment_arms) {
    arm_name <- if (multiarm) arm else "treated"
    
    # Subset to control + this arm
    idx <- W_factor %in% c(control, arm)
    W_binary <- as.integer(W_factor[idx] == arm)  # 1 = treatment arm, 0 = control
    X_sub <- X_df[idx, , drop = FALSE]
    X_matrix_sub <- X_matrix[idx, , drop = FALSE]
    Y_sub <- if (!is.null(Y)) Y[idx] else NULL
    
    n_per_arm[[arm_name]] <- list(
      control = sum(W_binary == 0),
      treated = sum(W_binary == 1)
    )
    
    # 1. Run fastcpt for formal balance test
    fastcpt_base_args <- list(
      Z = X_sub, 
      T = W_binary, 
      class.methods = class.method,
      perm.N = perm.N,
      alpha = alpha,
      R.seed = seed
    )
    fastcpt_all_args <- utils::modifyList(fastcpt_base_args, fastcpt.args)
    balance_test <- do.call(fastcpt, fastcpt_all_args)
    balance_tests[[arm_name]] <- balance_test
    passed_vec[arm_name] <- balance_test$pval > alpha
    
    # 2. Fit honest propensity models for visualization (boosted RF from grf)
    prop_real <- grf::boosted_regression_forest(
      X = X_matrix_sub, 
      Y = W_binary, 
      honesty = TRUE, 
      tune.parameters = "all", 
      seed = seed
    )
    pscores_real_list[[arm_name]] <- as.numeric(prop_real$predictions)
    
    # Null (permuted) treatment propensities
    set.seed(seed + which(treatment_arms == arm))
    W_null <- sample(W_binary)
    prop_null <- grf::boosted_regression_forest(
      X = X_matrix_sub, 
      Y = W_null, 
      honesty = TRUE, 
      tune.parameters = "all", 
      seed = seed
    )
    pscores_null_list[[arm_name]] <- as.numeric(prop_null$predictions)
    
    # 3. Optional treatment effect estimation (only if Y provided)
    if (!is.null(Y_sub)) {
      # Compute difference-in-means
      Y1 <- Y_sub[W_binary == 1]
      Y0 <- Y_sub[W_binary == 0]
      dim_est <- mean(Y1) - mean(Y0)
      dim_se <- sqrt(stats::var(Y1) / length(Y1) + stats::var(Y0) / length(Y0))
      dim_ci <- dim_est + c(-1, 1) * stats::qnorm(0.975) * dim_se
      
      dim_results[[arm_name]] <- list(
        estimate = dim_est,
        std.err = dim_se,
        ci = dim_ci
      )
      
      # Fit causal forest and get AIPW estimate (full model)
      cf <- grf::causal_forest(
        X = X_matrix_sub,
        Y = Y_sub,
        W = W_binary,
        seed = seed,
        num.trees = 2000
      )
      cf_list[[arm_name]] <- cf
      
      ate <- grf::average_treatment_effect(cf, target.sample = "all")
      
      aipw_results[[arm_name]] <- list(
        estimate = ate["estimate"],
        std.err = ate["std.err"],
        ci = ate["estimate"] + c(-1, 1) * stats::qnorm(0.975) * ate["std.err"]
      )
      
      # AIPW with constant propensity score (marginal probability)
      cf_const <- grf::causal_forest(
        X = X_matrix_sub,
        Y = Y_sub,
        W = W_binary,
        W.hat = rep(mean(W_binary), length(W_binary)),
        seed = seed,
        num.trees = 2000
      )
      
      ate_const <- grf::average_treatment_effect(cf_const, target.sample = "all")
      
      aipw_const_results[[arm_name]] <- list(
        estimate = ate_const["estimate"],
        std.err = ate_const["std.err"],
        ci = ate_const["estimate"] + c(-1, 1) * stats::qnorm(0.975) * ate_const["std.err"]
      )
    }
  }
  
  # Build result object
  # For binary treatment, unwrap single-element lists for backward compatibility
  if (!multiarm) {
    result <- list(
      balance_test = balance_tests[[1]],
      dim = if (length(dim_results) > 0) dim_results[[1]] else NULL,
      aipw = if (length(aipw_results) > 0) aipw_results[[1]] else NULL,
      aipw_const = if (length(aipw_const_results) > 0) aipw_const_results[[1]] else NULL,
      passed = unname(passed_vec[1]),
      alpha = alpha,
      cf = if (length(cf_list) > 0) cf_list[[1]] else NULL,
      pscores_real = pscores_real_list[[1]],
      pscores_null = pscores_null_list[[1]],
      n = length(W),
      n_treated = n_per_arm[[1]]$treated,
      n_control = n_per_arm[[1]]$control,
      control = control,
      arms = treatment_arms,
      multiarm = FALSE
    )
  } else {
    result <- list(
      balance_test = balance_tests,
      dim = if (length(dim_results) > 0) dim_results else NULL,
      aipw = if (length(aipw_results) > 0) aipw_results else NULL,
      aipw_const = if (length(aipw_const_results) > 0) aipw_const_results else NULL,
      passed = passed_vec,
      alpha = alpha,
      cf = if (length(cf_list) > 0) cf_list else NULL,
      pscores_real = pscores_real_list,
      pscores_null = pscores_null_list,
      n = length(W),
      n_per_arm = n_per_arm,
      control = control,
      arms = treatment_arms,
      multiarm = TRUE
    )
  }
  
  class(result) <- "balance"
  return(result)
}

#' @rdname balance
#' @export
print.balance <- function(x, ...) {
  cat("\nBalance Assessment\n")
  
  if (!x$multiarm) {
    # Binary treatment
    status <- ifelse(x$passed, "PASS", "FAIL")
    cat(sprintf("  Control: '%s'\n", x$control))
    cat(sprintf("  Balance test p-value: %.4f (%s)\n", x$balance_test$pval, status))
    if (!is.null(x$dim) && !is.null(x$aipw) && !is.null(x$aipw_const)) {
      cat(sprintf("  DiM estimate:         %.4f (SE: %.4f)\n", x$dim$estimate, x$dim$std.err))
      cat(sprintf("  Outcome adj. (CF):    %.4f (SE: %.4f)\n", x$aipw_const$estimate, x$aipw_const$std.err))
      cat(sprintf("  Doubly robust (CF):   %.4f (SE: %.4f)\n", x$aipw$estimate, x$aipw$std.err))
    } else {
      cat("  Outcome Y not provided: skipping treatment effect estimates.\n")
    }
  } else {
    # Multi-arm treatment
    cat(sprintf("  Multi-arm analysis (%d treatment arms vs control '%s')\n", 
                length(x$arms), x$control))
    cat("\n")
    
    for (arm in x$arms) {
      status <- ifelse(x$passed[arm], "PASS", "FAIL")
      cat(sprintf("  [%s vs %s]\n", arm, x$control))
      cat(sprintf("    Balance test p-value: %.4f (%s)\n", x$balance_test[[arm]]$pval, status))
      if (!is.null(x$dim) && !is.null(x$dim[[arm]])) {
        cat(sprintf("    DiM estimate:         %.4f (SE: %.4f)\n", 
                    x$dim[[arm]]$estimate, x$dim[[arm]]$std.err))
        cat(sprintf("    Outcome adj. (CF):    %.4f (SE: %.4f)\n", 
                    x$aipw_const[[arm]]$estimate, x$aipw_const[[arm]]$std.err))
        cat(sprintf("    Doubly robust (CF):   %.4f (SE: %.4f)\n", 
                    x$aipw[[arm]]$estimate, x$aipw[[arm]]$std.err))
      } else {
        cat("    Outcome Y not provided: skipping treatment effect estimates.\n")
      }
      cat("\n")
    }
  }
  
  cat("  Use summary() for details, plot() to visualize.\n\n")
  invisible(x)
}

#' @param object A balance result object (for summary method).
#' @rdname balance
#' @export
summary.balance <- function(object, ...) {
  
  cat("\n")
  cat("========================================================================\n")
  cat("                    COVARIATE BALANCE ASSESSMENT                        \n")
  cat("========================================================================\n\n")
  
  # Section 1: Sample
  cat("1. SAMPLE CHARACTERISTICS\n")
  cat("------------------------------------------------------------------------\n")
  cat(sprintf("   Total observations:           %d\n", object$n))
  cat(sprintf("   Control group ('%s'):       ", object$control))
  

  if (!object$multiarm) {
    # Binary treatment
    cat(sprintf("%d (%.1f%%)\n", object$n_control, 100 * object$n_control / object$n))
    cat(sprintf("   Treatment group:              %d (%.1f%%)\n", 
                object$n_treated, 100 * object$n_treated / object$n))
  } else {
    # Multi-arm treatment
    n_control <- object$n_per_arm[[1]]$control
    cat(sprintf("%d (%.1f%%)\n", n_control, 100 * n_control / object$n))
    for (arm in object$arms) {
      n_arm <- object$n_per_arm[[arm]]$treated
      cat(sprintf("   Treatment arm '%s':       %d (%.1f%%)\n", 
                  arm, n_arm, 100 * n_arm / object$n))
    }
  }
  cat("\n")
  
  # Helper function to print results for one arm
  print_arm_results <- function(arm_name, balance_test, ps_real, ps_null, 
                                 dim_res, aipw_res, aipw_const_res, cf, passed, 
                                 alpha, section_prefix = "") {
    status <- ifelse(passed, "PASS", "FAIL")
    
    # Balance Test
    cat(sprintf("%s2. CLASSIFICATION PERMUTATION TEST\n", section_prefix))
    cat("------------------------------------------------------------------------\n")
    cat(sprintf("   Classifier:                   %s\n", 
                paste(balance_test$class.methods, collapse = ", ")))
    cat(sprintf("   Number of permutations:       %d\n", balance_test$perm.N))
    cat(sprintf("   Test statistic (observed):    %.4f\n", balance_test$teststat))
    cat(sprintf("   Null distribution mean:       %.4f\n", mean(balance_test$nulldist)))
    cat(sprintf("   Null distribution SD:         %.4f\n", stats::sd(balance_test$nulldist)))
    cat(sprintf("   P-value:                      %.4f\n", balance_test$pval))
    cat(sprintf("   Significance level (alpha):   %.2f\n", alpha))
    cat(sprintf("   Result:                       %s\n", status))
    cat("\n")
    
    # Propensity Score Diagnostics
    cat(sprintf("%s3. PROPENSITY SCORE DIAGNOSTICS (Boosted Regression Forest)\n", section_prefix))
    cat("------------------------------------------------------------------------\n")
    cat("   Real treatment assignment:\n")
    cat(sprintf("      Mean:                      %.4f\n", mean(ps_real)))
    cat(sprintf("      SD:                        %.4f\n", stats::sd(ps_real)))
    cat(sprintf("      Range:                     [%.4f, %.4f]\n", min(ps_real), max(ps_real)))
    cat("   Null (permuted) assignment:\n")
    cat(sprintf("      Mean:                      %.4f\n", mean(ps_null)))
    cat(sprintf("      SD:                        %.4f\n", stats::sd(ps_null)))
    cat(sprintf("      Range:                     [%.4f, %.4f]\n", min(ps_null), max(ps_null)))
    cat("   Distributional comparison:\n")
    cat(sprintf("      Difference in means:       %.4f\n", mean(ps_real) - mean(ps_null)))
    cat(sprintf("      Ratio of SDs:              %.4f\n", stats::sd(ps_real) / stats::sd(ps_null)))
    cat("\n")
    
    # Treatment Effects
    cat(sprintf("%s4. TREATMENT EFFECT ESTIMATES\n", section_prefix))
    cat("------------------------------------------------------------------------\n")
    if (!is.null(dim_res) && !is.null(aipw_res) && !is.null(aipw_const_res)) {
      cat("   Difference-in-Means (unadjusted):\n")
      cat(sprintf("      Estimate:                  %.4f\n", dim_res$estimate))
      cat(sprintf("      Standard error:            %.4f\n", dim_res$std.err))
      cat(sprintf("      95%% CI (normal approx.):  [%.4f, %.4f]\n", 
                  dim_res$ci[1], dim_res$ci[2]))
      cat("\n")
      cat("   Outcome-adjusted (causal forest, no propensity weighting):\n")
      cat(sprintf("      Estimate:                  %.4f\n", aipw_const_res$estimate))
      cat(sprintf("      Standard error:            %.4f\n", aipw_const_res$std.err))
      cat(sprintf("      95%% CI (inf. jackknife):  [%.4f, %.4f]\n", 
                  aipw_const_res$ci[1], aipw_const_res$ci[2]))
      cat("\n")
      cat("   Doubly robust (causal forest with propensity weighting):\n")
      cat(sprintf("      Estimate:                  %.4f\n", aipw_res$estimate))
      cat(sprintf("      Standard error:            %.4f\n", aipw_res$std.err))
      cat(sprintf("      95%% CI (inf. jackknife):  [%.4f, %.4f]\n", 
                  aipw_res$ci[1], aipw_res$ci[2]))
      if (!is.null(cf)) {
        cat(sprintf("      Number of trees:           %d\n", cf$`_num_trees`))
      }
      cat("\n")
    } else {
      cat("   Outcome Y not provided: skipping treatment effect estimates.\n")
      cat("\n")
    }
  }
  
  if (!object$multiarm) {
    # Binary treatment - single output
    print_arm_results(
      arm_name = "treated",
      balance_test = object$balance_test,
      ps_real = object$pscores_real,
      ps_null = object$pscores_null,
      dim_res = object$dim,
      aipw_res = object$aipw,
      aipw_const_res = object$aipw_const,
      cf = object$cf,
      passed = object$passed,
      alpha = object$alpha
    )
    
    # Interpretation
    cat("5. INTERPRETATION\n")
    cat("------------------------------------------------------------------------\n")
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
    
  } else {
    # Multi-arm treatment - loop over arms
    for (i in seq_along(object$arms)) {
      arm <- object$arms[i]
      cat("========================================================================\n")
      cat(sprintf("  COMPARISON: %s vs %s (Control)\n", arm, object$control))
      cat("========================================================================\n\n")
      
      print_arm_results(
        arm_name = arm,
        balance_test = object$balance_test[[arm]],
        ps_real = object$pscores_real[[arm]],
        ps_null = object$pscores_null[[arm]],
        dim_res = if (!is.null(object$dim)) object$dim[[arm]] else NULL,
        aipw_res = if (!is.null(object$aipw)) object$aipw[[arm]] else NULL,
        aipw_const_res = if (!is.null(object$aipw_const)) object$aipw_const[[arm]] else NULL,
        cf = if (!is.null(object$cf)) object$cf[[arm]] else NULL,
        passed = object$passed[arm],
        alpha = object$alpha
      )
    }
    
    # Overall interpretation for multi-arm
    cat("========================================================================\n")
    cat("  OVERALL INTERPRETATION\n")
    cat("========================================================================\n\n")
    n_passed <- sum(object$passed)
    n_total <- length(object$passed)
    if (all(object$passed)) {
      cat(sprintf("   All %d pairwise balance tests PASSED.\n", n_total))
      cat("   The classifier cannot distinguish any treatment arm from control\n")
      cat("   better than chance.\n")
    } else if (!any(object$passed)) {
      cat(sprintf("   All %d pairwise balance tests FAILED.\n", n_total))
      cat("   The classifier can distinguish all treatment arms from control.\n")
    } else {
      cat(sprintf("   %d of %d pairwise balance tests passed.\n", n_passed, n_total))
      cat(sprintf("   Failed comparisons: %s\n", 
                  paste(names(object$passed)[!object$passed], collapse = ", ")))
    }
    cat("\n")
  }
  
  # Estimator interpretation (only if Y was provided)
  has_effects <- if (object$multiarm) !is.null(object$dim) else (!is.null(object$dim) && !is.null(object$aipw))
  if (has_effects) {
    cat("   ESTIMATOR GUIDE:\n")
    cat("   - DiM: Simple difference in means, no covariate adjustment.\n")
    cat("   - Outcome-adjusted: Uses causal forest for outcome regression only.\n")
    cat("   - Doubly robust: Combines outcome regression with propensity weighting.\n")
    cat("   If balance passes, all three should be similar. Large differences\n")
    cat("   may indicate model misspecification or residual confounding.\n")
    cat("\n")
  }
  
  invisible(object)
}

#' Plot method for balance objects
#'
#' @param x A balance result object.
#' @param which Character vector specifying which plots to create. Options are "pscores", "null_dist", "effects", or "all".
#' @param combined Logical. If TRUE, displays all three plots in a combined panel. Default is TRUE. 
#' @param breaks Number of breaks for histograms. Default is 15.
#' @param ... Additional arguments (currently unused).
#' @rdname balance
#' @export
plot.balance <- function(x, which = "all", combined = TRUE, breaks = 15, ...) {
  
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
  
  has_effects <- !is.null(x$dim) && !is.null(x$aipw) && !is.null(x$aipw_const)
  
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
          x = "Test Statistic (Classification Accuracy)",
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
      col_dim <- "#004225"       # British racing green
      col_aipw_const <- "#5ab0c0" # Teal (matching pdp)
      col_aipw <- "#FF8C00"      # Dark orange
      
      # Build data for ggdist (three estimators)
      effect_df <- data.frame(
        estimator = factor(c("Difference\nin Means", "Outcome\nAdjusted", "Doubly\nRobust"), 
                           levels = c("Difference\nin Means", "Outcome\nAdjusted", "Doubly\nRobust")),
        estimate = c(x$dim$estimate, x$aipw_const$estimate, x$aipw$estimate),
        se = c(x$dim$std.err, x$aipw_const$std.err, x$aipw$std.err)
      )
      
      # Create distributional data for ggdist
      effect_df$.dist <- lapply(seq_len(nrow(effect_df)), function(i) {
        distributional::dist_normal(effect_df$estimate[i], effect_df$se[i])
      })
      effect_df$.dist <- do.call(c, effect_df$.dist)
      
      plots$effects <- ggplot2::ggplot(effect_df, ggplot2::aes(x = estimator, ydist = .dist, color = estimator)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        ggdist::stat_pointinterval(
          .width = c(0.80, 0.95, 0.99),
          point_size = 3,
          interval_size_range = c(0.8, 2.5)
        ) +
        ggplot2::scale_color_manual(values = c(col_dim, col_aipw_const, col_aipw)) +
        g_theme() +
        ggplot2::labs(
          title = "C. Treatment Effect Estimates",
          x = NULL,
          y = "Effect Estimate",
          caption = "Note: Point estimates with 80%, 95%, 99% CIs. Intervals based on normal approximation."
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(size = 9, color = "black"),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank()
        )
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
    
    # Panel B: Classification permutation test null distribution (faceted by arm)
    if ("null_dist" %in% which) {
      # Build combined data frame for all arms
      null_df_list <- lapply(x$arms, function(arm) {
        data.frame(
          arm = sprintf("%s vs %s", arm, x$control),
          val = x$balance_test[[arm]]$nulldist
        )
      })
      null_df <- do.call(rbind, null_df_list)
      null_df$arm <- factor(null_df$arm, levels = unique(null_df$arm))
      
      # Test statistics and p-values per arm
      test_df <- do.call(rbind, lapply(x$arms, function(arm) {
        bt <- x$balance_test[[arm]]
        test_pass <- bt$pval > x$alpha
        data.frame(
          arm = sprintf("%s vs %s", arm, x$control),
          testval = bt$teststat,
          pval = bt$pval,
          passed = test_pass,
          status_text = ifelse(test_pass, "Pass", "Fail"),
          color = ifelse(test_pass, col_pass, col_fail)
        )
      }))
      test_df$arm <- factor(test_df$arm, levels = levels(null_df$arm))
      
      # Calculate axis limits across all arms
      all_nulldist <- unlist(lapply(x$arms, function(arm) x$balance_test[[arm]]$nulldist))
      all_testvals <- sapply(x$arms, function(arm) x$balance_test[[arm]]$teststat)
      x_min <- floor(min(all_nulldist) * 10) / 10
      x_max <- ceiling(max(c(all_nulldist, all_testvals)) * 10) / 10
      
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
        ggplot2::geom_vline(data = test_df, ggplot2::aes(xintercept = testval, color = color), 
                            linewidth = 1.5, show.legend = FALSE) +
        ggplot2::geom_text(data = test_df, 
                           ggplot2::aes(x = testval, y = 0.95, 
                                        label = sprintf("p = %.3f (%s)", pval, status_text),
                                        color = color),
                           hjust = -0.1, size = 3, fontface = "bold", show.legend = FALSE) +
        ggplot2::scale_color_identity() +
        ggplot2::facet_wrap(~ arm) +
        g_theme() +
        ggplot2::labs(
          title = "B. Classification Permutation Test by Treatment Arm",
          x = "Test Statistic (Classification Accuracy)",
          y = "Density",
          caption = sprintf("Note: Classifier = %s, alpha = %.2f.", 
                            paste(x$balance_test[[1]]$class.methods, collapse = ", "), x$alpha)
        ) +
        ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0))
    }
    
    # Panel C: Treatment effect estimates (faceted by arm)
    if ("effects" %in% which && has_effects) {
      # Colors for estimators
      col_dim <- "#004225"       # British racing green
      col_aipw_const <- "#5ab0c0" # Teal
      col_aipw <- "#FF8C00"      # Dark orange
      
      # Build combined data frame for all arms
      effect_df_list <- lapply(x$arms, function(arm) {
        data.frame(
          arm = sprintf("%s vs %s", arm, x$control),
          estimator = factor(c("Difference\nin Means", "Outcome\nAdjusted", "Doubly\nRobust"), 
                             levels = c("Difference\nin Means", "Outcome\nAdjusted", "Doubly\nRobust")),
          estimate = c(x$dim[[arm]]$estimate, x$aipw_const[[arm]]$estimate, x$aipw[[arm]]$estimate),
          se = c(x$dim[[arm]]$std.err, x$aipw_const[[arm]]$std.err, x$aipw[[arm]]$std.err)
        )
      })
      effect_df <- do.call(rbind, effect_df_list)
      effect_df$arm <- factor(effect_df$arm, levels = unique(effect_df$arm))
      
      # Create distributional data for ggdist
      effect_df$.dist <- lapply(seq_len(nrow(effect_df)), function(i) {
        distributional::dist_normal(effect_df$estimate[i], effect_df$se[i])
      })
      effect_df$.dist <- do.call(c, effect_df$.dist)
      
      plots$effects <- ggplot2::ggplot(effect_df, ggplot2::aes(x = estimator, ydist = .dist, color = estimator)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        ggdist::stat_pointinterval(
          .width = c(0.80, 0.95, 0.99),
          point_size = 3,
          interval_size_range = c(0.8, 2.5)
        ) +
        ggplot2::scale_color_manual(values = c(col_dim, col_aipw_const, col_aipw)) +
        ggplot2::facet_wrap(~ arm) +
        g_theme() +
        ggplot2::labs(
          title = "C. Treatment Effect Estimates by Arm",
          x = NULL,
          y = "Effect Estimate",
          caption = "Note: Point estimates with 80%, 95%, 99% CIs. Intervals based on normal approximation."
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(size = 8, color = "black", angle = 45, hjust = 1),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank()
        )
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

