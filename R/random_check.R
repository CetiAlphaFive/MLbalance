# MAIN RANDOM CHECK FUNCTION FOR BINARY TREATMENTS ####
#
# Suppress R CMD check notes for non-standard evaluation in ggplot
utils::globalVariables(c("val", "var"))
#
# Check dependencies
if (!requireNamespace("grf", quietly = TRUE)) {
  stop(
    "Package \"grf\" must be installed to use MLbalance.",
    call. = FALSE
  )
}
#
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package \"ggplot2\" must be installed to use MLbalance.",
    call. = FALSE
  )
}
#
if (!requireNamespace("ggdist", quietly = TRUE)) {
  stop(
    "Package \"ggdist\" must be installed to use MLbalance.",
    call. = FALSE
  )
}
#
# This measure of variable importance is explained in the appendix, comes from grf. Function to arrange scores
#
#' Variable Importance Function
#' @import grf
#' @import ggdist
#' @import ggplot2
#' @param model Trained GRF Model Object
#' @description
#' This convenience function takes a trained grf model object and returns a data frame of variable importance scores using grf's simple variable importance metric.
#'
#' @examples x <- data.frame(X1 = rnorm(1000))
#' @examples y <- rnorm(1000)
#' @examples model <- grf::regression_forest(X = x,Y = y)
#' @examples vip(model)
#' @export
vip <- function(model){
  vip_scores <- data.frame(varname = colnames(model$X.orig),vip = grf::variable_importance(model))
  vip_scores[order(vip_scores$vip, decreasing = T),]
}
#
#' Balance Permutation Test
#
#' @param W_real Real treatment assignment vector.
#' @param W_sim Simulated treatment assignment vector. If not provided, permuted W_real is used.
#' @param X Pre-treatment covariate matrix or data frame.
#' @param R.seed Random seed used in set.seed (for replicability).
#' @param grf.seed Random seed used in grf's seed (for replicability).
#' @param breaks number of breaks in output histogram. Default is 15.
#' @param facet facet by treatment assignment. Default is FALSE.
#' @description
#' This is the main balance permutation test function. First, the function attempts to model treatment assignment (W_real) as a function of pre-treatment covariates (X). It does so using an honest, boosted random forest (see Ghosal and Hooker 2018) with built-in hyperparameter tuning. This model is used to generate real treatment propensity scores. Then, we build a second boosted random forest model using the same pre-treatment covariates and tuning parameter settings but with either simulated, or randomly permuted, treatment assignment as the outcome variable. The function proceeds to output both the real and null treatment propensity scores as well as diagnostics and a plot comparing the distributions.
#'
#' The purpose of this exercise is to compare the real treatment propensity distribution to a null distribution where treatment assignment is correctly orthogonal to pre-treatment covariates. To interpret the results, it's advisable to notice any extreme, deterministic treatment propensity scores near zero or one, or any other divergences from design expectations. In general, if randomization succeeded the two distributions should closely overlap with similar means and variances. If the results are at all unclear, it's advisable to estimate average treatment effects via a method that accounts for propensity to treatment (e.g., augmented inverse propensity weighting, overlap weighting, etc.).
#'
#'
#' @examples
#' \dontrun{
#' n <- 1000
#' p <- 20
#' X <- matrix(rnorm(n*p,0,1),n,p)
#' w_real <- rbinom(n, 1, ifelse(.021 + abs(.4*X[,4] - .5*X[,8]) < 1,
#'                  .021 + abs(.4*X[,4] - .5*X[,8]), 1))
#' df <- data.frame(w_real,X)
#' r.check <- random_check(W_real = df$w_real,
#'                         X = subset(df, select = -w_real))
#' r.check$plot
#' }
#' @export
random_check <- function(W_real, W_sim = NULL, X,R.seed = 1995, grf.seed = 1995, breaks = 15,facet = FALSE){

  #Set the seed, default is 1995
  set.seed(R.seed)

  #Print message if permutation selected
  if(is.null(W_sim) == TRUE){
    cat("No Simulated Assignment Vector Provided, Null Distribution Generated Using Permutated Treatment Assignment.\n\n\n")} else {
      cat("Simulated Assignment Vector Provided, Null Distribution Generated Using Simulated Treatment Assignment.\n\n\n")
    }

  #Print simple count table(s)
  if(length(unique(W_real)) <= 2 & !is.null(W_sim)){cat("Simple Count Table(s)\n\n"); print(table(W_real)); print(table(W_sim))}

  #Check inputs for correct formats
  if(is.vector(W_real) != TRUE)
    stop("Error: W_real must be a vector")

  if(is.vector(W_sim) != TRUE & is.null(W_sim) == FALSE)
    stop("Error: W_sim must be a vector or NULL.")

  if(is.matrix(X) != TRUE & is.data.frame(X) != TRUE)
    stop("Error: X must be a matrix or data frame.")

  #Check if simulated treatment assignments provided, if not, permute the real treatment assignment vector.
  if(is.null(W_sim) == TRUE){W_sim <- sample(W_real,size = length(W_real),replace = FALSE)}

  # Build a treatment propensity model with the real treatment assignment vector. Boosted reg forest from grf.
  g.real  <- grf::boosted_regression_forest(X = X, Y = W_real, honesty = T, tune.parameters = "all", seed = grf.seed)

  # Build a treatment propensity model with the simulated treatment assignment vector. Boosted reg forest from grf.
  g.sim   <- grf::boosted_regression_forest(X = X, Y = W_sim,  honesty = T, tune.parameters = "all", seed = grf.seed)

  # Build a data frame for the diagnostics plot
  plot.df <- data.frame(var = factor(c(rep("Real",NROW(g.real$predictions)),rep("Null",NROW(g.sim$predictions)))),
                        treat = c(W_real, W_sim),
                        val = c(g.real$predictions,g.sim$predictions))

  #ggplot theme, no package required
  g_theme <- function(){
    theme(plot.title = element_text(size=14, face="bold", hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "white", linewidth = 0.5, linetype = "solid"),
          axis.line = element_line(linewidth = .5, color = "black"),
          axis.title.x = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
          axis.text.x = element_text(color = "grey30", size=10),
          axis.text.y = element_text(color = "grey30", size=10),
          panel.grid.major.x = element_line(colour = "grey80"),
          plot.caption = element_text(hjust = 0),
          text=element_text(size=12,  family="serif"),
          legend.position = c(.1,.75),
          axis.ticks = element_blank(),
          legend.background=element_blank(),
          panel.spacing = unit(2, "lines"))
  }

  #Create the overlapping histogram ggplot2
  g <- ggplot(plot.df, aes(x = val,fill = var)) +
    ggdist::stat_histinterval(slab_color = "gray70",
                              outline_bars = TRUE,
                              alpha = .75,
                              point_alpha = 0,
                              slab_linewidth = .5,
                              breaks = breaks,
                              interval_alpha = 0) +
    geom_vline(xintercept = mean(g.real$predictions),color = "darkorange1",linetype = "dotdash",linewidth = .5) +
    geom_vline(xintercept = mean(g.sim$predictions),color = "dodgerblue1",linetype = "dotdash",linewidth = .5) +
    g_theme() +
    xlab("Treatment Propensity Scores") +
    ylab("Density") +
    labs(caption = expression(italic("Note: Dotted lines represent mean values for the null and real treatment propensity distributions."))) +
    #ggtitle("Randomization Check: Real Treatment Propensities vs. Permuted Null Distribution") +
    scale_fill_manual(values = c("dodgerblue1","darkorange1")) +
    guides(fill=guide_legend(title="")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    if(max(plot.df$val) > 1 | min(plot.df$val) < 0){scale_x_continuous(expand = c(0, 0))}else{scale_x_continuous(limits = c(0, 1.01), expand = c(0, 0))}

  g <- g +
    if(facet == TRUE){facet_wrap(~treat)}

  #Create the results object list
  results <- list(
    "prop.model.real" = g.real,
    "prop.model.real.tuning" = g.real$forests[[1]]$tunable.params,
    "treat.props.real" = g.real$predictions,
    "imp.predictors" = vip(g.real$forests[[1]]),
    "prop.model.sim" = g.sim,
    "prop.model.real.tuning" = g.sim$forests[[1]]$tunable.params,
    "treat.props.sim" = g.sim$predictions,
    "plot.df" = plot.df,
    "plot" = g)

  #return the results
  return(results)
}
#
