# Suppress R CMD check notes for non-standard evaluation
utils::globalVariables(c("estimator", "estimate", "var", "val", ".dist"))

#' Complete Balance Assessment and Treatment Effect Estimation
#'
#' @description
#' A unified function for covariate balance assessment and treatment effect estimation.
#' Combines a formal balance test (via classification permutation test), visual diagnostics
#' (propensity score distributions), and treatment effect estimates using both difference-in-means
#' and augmented inverse propensity weighting (AIPW).
#'
#' @param Y Outcome vector (numeric) or \code{NULL}. If \code{NULL}, treatment effect
#'   estimation (and the treatment effect plot) is skipped.
#' @param W Treatment assignment vector (binary: 0/1 or logical).
#' @param X Pre-treatment covariate matrix or data frame.
#' @param alpha Significance level for balance test. Default is 0.05.
#' @param perm.N Number of permutations for the balance test. Default is 1000.
#' @param class.method Classification method for balance test. Default is "ferns".
#' @param seed Random seed for reproducibility. Default is 1995.
#' @param fastcpt.args A named list of additional arguments to pass to \code{\link{fastcpt}}. For example, \code{fastcpt.args = list(parallel = TRUE, leaveout = 0.2)}.
#'
#' @return A list of class "balance" containing:
#' \item{balance_test}{Results from fastcpt including p-value and propensity scores.}
#' \item{dim}{Difference-in-means estimate with standard error and confidence interval (only if \code{Y} is provided).}
#' \item{aipw}{Doubly robust estimate from causal forest (with propensity weighting) with standard error and CI (only if \code{Y} is provided).}
#' \item{aipw_const}{Outcome-adjusted estimate from causal forest (no propensity weighting) with SE and CI (only if \code{Y} is provided).}
#' \item{passed}{Logical indicating whether the balance test passed.}
#' \item{alpha}{The significance level used.}
#' \item{cf}{The fitted causal_forest object for advanced users (only if \code{Y} is provided).}
#'
#' @examples
#' \dontrun{
#' # Generate example data
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
#' }
#'
#' @export
balance <- function(Y = NULL, W, X, alpha = 0.05, perm.N = 1000, class.method = "ferns", seed = 1995, 
                    fastcpt.args = list()) {
  
  # Check for grf
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required for propensity score modeling (and optional treatment effect estimation).", call. = FALSE)
  }
  
  # Input validation
  if (length(unique(W)) != 2) stop("W must be binary (two unique values).", call. = FALSE)
  if (length(W) != nrow(X)) {
    stop("W and X must have the same number of observations.", call. = FALSE)
  }
  if (!is.null(Y)) {
    if (!is.numeric(Y)) stop("Y must be a numeric vector.", call. = FALSE)
    if (length(Y) != length(W)) {
      stop("Y and W must have the same number of observations.", call. = FALSE)
    }
  }
  
  # Convert W to 0/1 if needed
  W <- as.integer(as.factor(W)) - 1L
  
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
  
  # 1. Run fastcpt for formal balance test
  fastcpt_base_args <- list(
    Z = X_df, 
    T = W, 
    class.methods = class.method,
    perm.N = perm.N,
    alpha = alpha,
    R.seed = seed
  )
  # Merge with user-supplied args (user args override defaults)
  fastcpt_all_args <- utils::modifyList(fastcpt_base_args, fastcpt.args)
  balance_test <- do.call(fastcpt, fastcpt_all_args)
  
  passed <- balance_test$pval > alpha
  
  # 2. Fit honest propensity models for visualization (boosted RF from grf)
  # Real treatment propensities
  prop_real <- grf::boosted_regression_forest(
    X = X_matrix, 
    Y = W, 
    honesty = TRUE, 
    tune.parameters = "all", 
    seed = seed
  )
  pscores_real <- prop_real$predictions
  
  # Null (permuted) treatment propensities
  # Reset seed for reproducibility (offset from main seed to ensure distinct from fastcpt permutations)
  set.seed(seed + 1L)
  W_null <- sample(W)
  prop_null <- grf::boosted_regression_forest(
    X = X_matrix, 
    Y = W_null, 
    honesty = TRUE, 
    tune.parameters = "all", 
    seed = seed
  )
  pscores_null <- prop_null$predictions
  
  # 3. Optional treatment effect estimation (only if Y provided)
  if (!is.null(Y)) {
    # Compute difference-in-means
    Y1 <- Y[W == 1]
    Y0 <- Y[W == 0]
    dim_est <- mean(Y1) - mean(Y0)
    dim_se <- sqrt(stats::var(Y1) / length(Y1) + stats::var(Y0) / length(Y0))
    dim_ci <- dim_est + c(-1, 1) * stats::qnorm(0.975) * dim_se
    
    dim_result <- list(
      estimate = dim_est,
      std.err = dim_se,
      ci = dim_ci
    )
    
    # Fit causal forest and get AIPW estimate (full model)
    cf <- grf::causal_forest(
      X = X_matrix,
      Y = Y,
      W = W,
      seed = seed,
      num.trees = 2000
    )
    
    ate <- grf::average_treatment_effect(cf, target.sample = "all")
    
    aipw_result <- list(
      estimate = ate["estimate"],
      std.err = ate["std.err"],
      ci = ate["estimate"] + c(-1, 1) * stats::qnorm(0.975) * ate["std.err"]
    )
    
    # AIPW with constant propensity score (marginal probability)
    cf_const <- grf::causal_forest(
      X = X_matrix,
      Y = Y,
      W = W,
      W.hat = rep(mean(W), length(W)),  # Constant propensity = marginal prob
      seed = seed,
      num.trees = 2000
    )
    
    ate_const <- grf::average_treatment_effect(cf_const, target.sample = "all")
    
    aipw_const_result <- list(
      estimate = ate_const["estimate"],
      std.err = ate_const["std.err"],
      ci = ate_const["estimate"] + c(-1, 1) * stats::qnorm(0.975) * ate_const["std.err"]
    )
  } else {
    dim_result <- NULL
    aipw_result <- NULL
    aipw_const_result <- NULL
    cf <- NULL
  }
  
  # Build result object
  result <- list(
    balance_test = balance_test,
    dim = dim_result,
    aipw = aipw_result,
    aipw_const = aipw_const_result,
    passed = passed,
    alpha = alpha,
    cf = cf,
    pscores_real = as.numeric(pscores_real),
    pscores_null = as.numeric(pscores_null),
    n = length(W),
    n_treated = sum(W),
    n_control = sum(W == 0)
  )
  
  class(result) <- "balance"
  return(result)
}

#' @rdname balance
#' @export
print.balance <- function(x, ...) {
  status <- ifelse(x$passed, "PASS", "FAIL")
  
  cat("\nBalance Assessment\n")
  cat(sprintf("  Balance test p-value: %.4f (%s)\n", x$balance_test$pval, status))
  if (!is.null(x$dim) && !is.null(x$aipw) && !is.null(x$aipw_const)) {
    cat(sprintf("  DiM estimate:         %.4f (SE: %.4f)\n", x$dim$estimate, x$dim$std.err))
    cat(sprintf("  Outcome adj. (CF):    %.4f (SE: %.4f)\n", x$aipw_const$estimate, x$aipw_const$std.err))
    cat(sprintf("  Doubly robust (CF):   %.4f (SE: %.4f)\n", x$aipw$estimate, x$aipw$std.err))
  } else {
    cat("  Outcome Y not provided: skipping treatment effect estimates.\n")
  }
  cat("  Use summary() for details, plot() to visualize.\n\n")
  
  invisible(x)
}

#' @param object A balance result object (for summary method).
#' @rdname balance
#' @export
summary.balance <- function(object, ...) {
  status <- ifelse(object$passed, "PASS", "FAIL")
  
  # Propensity score diagnostics
  ps_real <- object$pscores_real
  ps_null <- object$pscores_null
  
  cat("\n")
  cat("========================================================================\n")
  cat("                    COVARIATE BALANCE ASSESSMENT                        \n")
  cat("========================================================================\n\n")
  
  # Section 1: Sample
  cat("1. SAMPLE CHARACTERISTICS\n")
  cat("------------------------------------------------------------------------\n")
  cat(sprintf("   Total observations:           %d\n", object$n))
  cat(sprintf("   Treatment group (W=1):        %d (%.1f%%)\n", 
              object$n_treated, 100 * object$n_treated / object$n))
  cat(sprintf("   Control group (W=0):          %d (%.1f%%)\n", 
              object$n_control, 100 * object$n_control / object$n))
  cat("\n")
  
  # Section 2: Balance Test
  cat("2. CLASSIFICATION PERMUTATION TEST\n")
  cat("------------------------------------------------------------------------\n")
  cat(sprintf("   Classifier:                   %s\n", 
              paste(object$balance_test$class.methods, collapse = ", ")))
  cat(sprintf("   Number of permutations:       %d\n", object$balance_test$perm.N))
  cat(sprintf("   Test statistic (observed):    %.4f\n", object$balance_test$teststat))
  cat(sprintf("   Null distribution mean:       %.4f\n", mean(object$balance_test$nulldist)))
  cat(sprintf("   Null distribution SD:         %.4f\n", stats::sd(object$balance_test$nulldist)))
  cat(sprintf("   P-value:                      %.4f\n", object$balance_test$pval))
  cat(sprintf("   Significance level (alpha):   %.2f\n", object$alpha))
  cat(sprintf("   Result:                       %s\n", status))
  cat("\n")
  
  # Section 3: Propensity Score Diagnostics
  cat("3. PROPENSITY SCORE DIAGNOSTICS (Boosted Regression Forest)\n")
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
  
  # Section 4: Treatment Effects (optional)
  cat("4. TREATMENT EFFECT ESTIMATES\n")
  cat("------------------------------------------------------------------------\n")
  if (!is.null(object$dim) && !is.null(object$aipw) && !is.null(object$aipw_const)) {
    cat("   Difference-in-Means (unadjusted):\n")
    cat(sprintf("      Estimate:                  %.4f\n", object$dim$estimate))
    cat(sprintf("      Standard error:            %.4f\n", object$dim$std.err))
    cat(sprintf("      95%% CI (normal approx.):  [%.4f, %.4f]\n", 
                object$dim$ci[1], object$dim$ci[2]))
    cat("\n")
    cat("   Outcome-adjusted (causal forest, no propensity weighting):\n")
    cat(sprintf("      Estimate:                  %.4f\n", object$aipw_const$estimate))
    cat(sprintf("      Standard error:            %.4f\n", object$aipw_const$std.err))
    cat(sprintf("      95%% CI (inf. jackknife):  [%.4f, %.4f]\n", 
                object$aipw_const$ci[1], object$aipw_const$ci[2]))
    cat("\n")
    cat("   Doubly robust (causal forest with propensity weighting):\n")
    cat(sprintf("      Estimate:                  %.4f\n", object$aipw$estimate))
    cat(sprintf("      Standard error:            %.4f\n", object$aipw$std.err))
    cat(sprintf("      95%% CI (inf. jackknife):  [%.4f, %.4f]\n", 
                object$aipw$ci[1], object$aipw$ci[2]))
    if (!is.null(object$cf)) {
      cat(sprintf("      Number of trees:           %d\n", object$cf$`_num_trees`))
    }
    cat("\n")
  } else {
    cat("   Outcome Y not provided: skipping treatment effect estimates.\n")
    cat("\n")
  }
  
  # Section 5: Interpretation
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
  
  # Estimator interpretation (only if Y was provided)
  if (!is.null(object$dim) && !is.null(object$aipw) && !is.null(object$aipw_const)) {
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
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      panel.spacing = ggplot2::unit(2, "lines"),
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
  
  # Combined multi-panel display
  if (combined && length(plots) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required for combined plots.", call. = FALSE)
    }
    return(patchwork::wrap_plots(plotlist = unname(plots)))
  }
  
  if (length(plots) == 1) {
    return(plots[[1]])
  }
  
  return(plots)
}

