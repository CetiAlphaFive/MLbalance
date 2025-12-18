# Suppress R CMD check notes for non-standard evaluation in ggplot
utils::globalVariables(c("null"))

#' Plot fastcpt Results
#'
#' @description
#' Creates a visualization of the Classification Permutation Test results. Displays the permuted null distribution as a histogram with the observed test statistic overlaid as a vertical line. The plot indicates whether the test passed (p > 0.05) or failed.
#'
#' @param x A fastcpt result object as returned by \code{\link{fastcpt}}.
#' @param breaks Number of breaks for the histogram. Default is 25.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object displaying the null distribution histogram with the observed test statistic.
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
#'
#' # Plot the results
#' plot(result)
#' }
#'
#' @method plot fastcpt
#' @export
plot.fastcpt <- function(x, breaks = 25, ...){
  # package checks
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package \"ggdist\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # extract data and setup colors for figure
  alpha <- if (!is.null(x$alpha)) x$alpha else 0.05
  test_pass <- x$pval > alpha
  color_select <- ifelse(test_pass, "green", "red")
  subtitle_text <- ifelse(test_pass, "Pass", "Fail")
  testval <- x$teststat
  df <- data.frame(null = x$nulldist)

  # dynamic x-axis limits (rounded to nearest 0.1)
  x_min <- floor(min(x$nulldist) * 10) / 10
  x_max <- ceiling(max(c(x$nulldist, testval)) * 10) / 10

  # same theme as the random_check function
  g_theme <- function(){
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5, color = color_select),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white", linewidth = 0.5, linetype = "solid"),
      axis.line = ggplot2::element_line(linewidth = .5, color = "black"),
      axis.title.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(color = "grey30", size = 10),
      axis.text.y = ggplot2::element_text(color = "grey30", size = 10),
      panel.grid.major.x = ggplot2::element_line(colour = "grey80"),
      plot.caption = ggplot2::element_text(hjust = 0),
      text = ggplot2::element_text(size = 12, family = "serif"),
      legend.position = c(.1, .75),
      axis.ticks = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(2, "lines")
    )
  }

  # the Gagnon plot
  ggplot2::ggplot(df, ggplot2::aes(x = null)) +
    ggdist::stat_histinterval(
      slab_color = "gray70",
      outline_bars = TRUE,
      alpha = .75,
      point_alpha = 0,
      slab_linewidth = .5,
      breaks = breaks,
      interval_alpha = 0
    ) +
    ggplot2::geom_vline(xintercept = testval, color = color_select, linewidth = 1.25) +
    ggplot2::labs(
      x = "Permuted Null Distribution \n (Treatment Classification Model Fit)",
      y = "Density",
      title = "Classification Permutation Test Result",
      subtitle = subtitle_text
    ) +
    g_theme() +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
}

#' Summary of fastcpt Results
#'
#' @description
#' Prints a formatted summary of the Classification Permutation Test results including p-values, test statistics, and pass/fail status.
#'
#' @param object A fastcpt result object as returned by \code{\link{fastcpt}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the original object.
#'
#' @examples
#' \dontrun{
#' result <- fastcpt(Z, T, class.methods = "forest", perm.N = 100)
#' summary(result)
#' }
#'
#' @method summary fastcpt
#' @export
summary.fastcpt <- function(object, ...) {


  # Classifier name mapping
  classifier_names <- c(
    "ferns" = "Random Ferns",
    "forest" = "Random Forest",
    "lda" = "Linear Discriminant Analysis",
    "logistic" = "Logistic Regression",
    "logistic2" = "Logistic Regression (2-way)",
    "logistic2fast" = "Fast Logistic (2-way)",
    "glmnet" = "Elastic Net (glmnet)",
    "glmnet2" = "Elastic Net (2-way)",
    "ensemble" = "Ensemble"
  )

  # Metric name mapping
  metric_names <- c(
    "probability" = "Probability",
    "rate" = "Classification Rate",
    "mse" = "Mean Squared Error",
    "logscore" = "Log Score"
  )

  # Get full names
  get_classifier_name <- function(x) {
    if (x %in% names(classifier_names)) classifier_names[x] else x
  }

  get_metric_name <- function(x) {
    if (is.character(x) && x %in% names(metric_names)) metric_names[x] else "Probability"
  }

  cat("\n")
  cat("======================================================\n")
  cat("       Classification Permutation Test Results        \n")
  cat("======================================================\n\n")

  # Determine pass/fail
  alpha <- if (!is.null(object$alpha)) object$alpha else 0.05
  pass <- object$pval > alpha
  status <- ifelse(pass, "PASS", "FAIL")

  # Number of permutations
  n_perm <- if (!is.null(object$perm.N)) object$perm.N else length(object$nulldist)

  # Metric used
  metric_display <- get_metric_name(object$metric)

  # Single classifier case
  if (length(object$pvals) == 1) {
    classifier_display <- get_classifier_name(object$class.methods[1])
    cat(sprintf("  Classifier:        %s\n", classifier_display))
    cat(sprintf("  Metric:            %s\n", metric_display))
    cat(sprintf("  Test Statistic:    %.4f\n", object$teststat))
    cat(sprintf("  P-value:           %.4f\n", object$pval))
    cat(sprintf("  Permutations:      %d\n", n_perm))
    cat(sprintf("  Result:            %s (alpha = %.2f)\n", status, alpha))
  } else {
    # Multiple classifiers
    cat(sprintf("  Metric:            %s\n", metric_display))
    cat(sprintf("  Combined P-value:  %.4f\n", object$pval))
    cat(sprintf("  Permutations:      %d\n", n_perm))
    cat(sprintf("  Result:            %s (alpha = %.2f)\n", status, alpha))
    cat("\n")
    cat("  Individual Classifiers:\n")
    cat("  -------------------------------------------------\n")

    # Build results table
    methods <- names(object$pvals)
    for (i in seq_along(methods)) {
      if (methods[i] != "ensemble") {
        full_name <- get_classifier_name(methods[i])
        cat(sprintf("  %-28s  stat: %.4f  p: %.4f\n",
            full_name, object$teststat[i], object$pvals[i]))
      }
    }
    if ("ensemble" %in% methods) {
      cat(sprintf("  %-28s  stat: %.4f  p: %.4f\n",
          "Ensemble", object$teststat[length(object$teststat)],
          object$pvals["ensemble"]))
    }
  }

  cat("\n======================================================\n")
  if (pass) {
    cat("  Interpretation: Treatment groups appear BALANCED.\n")
    cat("  The classifier cannot distinguish between groups\n")
    cat(sprintf("  better than chance (p > %.2f).\n", alpha))
  } else {
    cat("  Interpretation: Treatment groups appear IMBALANCED.\n")
    cat("  The classifier can distinguish between groups\n")
    cat(sprintf("  better than chance (p <= %.2f).\n", alpha))
  }
  cat("======================================================\n\n")

  invisible(object)
}

#' Print fastcpt Results
#'
#' @description
#' Prints a brief overview of the Classification Permutation Test results.
#'
#' @param x A fastcpt result object as returned by \code{\link{fastcpt}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the original object.
#'
#' @method print fastcpt
#' @export
print.fastcpt <- function(x, ...) {
  alpha <- if (!is.null(x$alpha)) x$alpha else 0.05
  pass <- x$pval > alpha
  status <- ifelse(pass, "PASS", "FAIL")

  cat("\nClassification Permutation Test\n")
  cat(sprintf("  P-value: %.4f (%s)\n", x$pval, status))
  cat("  Use summary() for details, plot() to visualize.\n\n")

  invisible(x)
}
