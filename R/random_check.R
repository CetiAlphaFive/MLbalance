# MAIN RANDOM CHECK FUNCTION FOR BINARY TREATMENTS ####
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
#' @examples vip(grf_model_object)
#' @export
vip <- function(model){
  vip_scores <- data.frame(varname = colnames(model$X.orig),vip = grf::variable_importance(model))
  vip_scores[order(vip_scores$vip, decreasing = T),]
}
#
#' Permutation Balance Test
#
#' @param W_real Real Treatment Assignment Vector
#' @param W_sim Simulated Treatment Assignment Vector
#' @param R.seed Random seed used in set.seed
#' @param grf.seed Random seed used in grf's seed
#' @param breaks number of breaks in output histogram
#' @param facet facet by treatment assignment,default is FALSE
#' @examples n <- 1000
#' @examples p <- 20
#' @examples X <- matrix(rnorm(n*p,0,1),n,p)
#' @examples w_real <- rbinom(n, 1, ifelse(.021 + abs(.4*X[,4] - .5*X[,8]) < 1, .021 + abs(.4*X[,4] - .5*X[,8]), 1))
#' @examples df <- data.frame(w_real,X)
#' @examples r.check <- random_check(W_real = df$w_real, W_sim  = df$w_sim,X = subset(df,select = -w_real)); r.check
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

  # Build a treatment propensity model with the simulated treatment assignment vector. Lock tuning parameters to real model.
  g.sim   <- grf::boosted_regression_forest(X = X, Y = W_sim, honesty = T, seed = grf.seed,
                                            sample.fraction      = g.real$forests[[1]]$tunable.params$sample.fraction,
                                            mtry                 = g.real$forests[[1]]$tunable.params$mtry,
                                            min.node.size        = g.real$forests[[1]]$tunable.params$min.node.size,
                                            honesty.fraction     = g.real$forests[[1]]$tunable.params$honesty.fraction,
                                            honesty.prune.leaves = g.real$forests[[1]]$tunable.params$honesty.prune.leaves,
                                            alpha                = g.real$forests[[1]]$tunable.params$alpha,
                                            boost.steps          = length(g.real$forests))

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
    "treat.props" = g.real$predictions,
    "imp.predictors" = vip(g.real$forests[[1]]),
    "prop.model.sim" = g.sim,
    "prop.model.real.tuning" = g.sim$forests[[1]]$tunable.params,
    "treatment.props.sim" = g.sim$predictions,
    "plot.df" = plot.df,
    "plot" = g)

  # Add marker for detection of extreme treatment propensities
  results$extreme.props <- if(any(results$treat.props > .9 | results$treat.props < .1)){
    cat("\n\nWarning: Extreme Propensity scores detected (greater than .9 or less than .1).
                                      Examine $treat.props for more detail. \n\n")} else {
                                        cat("\n\nNo Extreme Propensity scores detected (greater than .9 or less than .1). \n\n")
                                      }

  #Print the plot by default
  results$g

  #return the results
  return(results)
}
#

