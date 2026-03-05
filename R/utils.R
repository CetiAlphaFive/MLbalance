# Internal helpers shared across plot functions

# Suppress R CMD check notes
utils::globalVariables(character(0))

#' @keywords internal
#' @noRd
.save_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv())
  else NULL
}

#' @keywords internal
#' @noRd
.restore_rng_state <- function(old_seed) {
  if (is.null(old_seed))
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  else
    assign(".Random.seed", old_seed, envir = globalenv())
}

#' @keywords internal
#' @noRd
.g_theme <- function(subtitle_color = NULL) {
  base <- ggplot2::theme(
    plot.title        = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    panel.background  = ggplot2::element_rect(fill = "white", colour = "white",
                                              linewidth = 0.5, linetype = "solid"),
    axis.line         = ggplot2::element_line(linewidth = 0.5, color = "black"),
    axis.title.x      = ggplot2::element_text(size = 12, face = "bold"),
    axis.title.y      = ggplot2::element_text(size = 12, face = "bold"),
    axis.text.x       = ggplot2::element_text(color = "grey30", size = 10),
    axis.text.y       = ggplot2::element_text(color = "grey30", size = 10),
    panel.grid.major.x = ggplot2::element_line(colour = "grey80"),
    plot.caption      = ggplot2::element_text(hjust = 0),
    text              = ggplot2::element_text(size = 12, family = "serif"),
    legend.position   = c(0.1, 0.75),
    axis.ticks        = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    panel.spacing     = ggplot2::unit(2, "lines")
  )

  if (!is.null(subtitle_color)) {
    base <- base + ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 12, face = "bold",
                                            hjust = 0.5, color = subtitle_color)
    )
  }

  base
}

#' @keywords internal
#' @noRd
.make_ci <- function(estimate, se, alpha = 0.05) {
  estimate + c(-1, 1) * stats::qnorm(1 - alpha / 2) * se
}

#' @keywords internal
#' @noRd
.validate_clusters_blocks <- function(clusters, blocks, n, paired) {
  if (!is.null(clusters)) {
    if (length(clusters) != n)
      stop("'clusters' must have the same length as the number of observations.", call. = FALSE)
    if (any(is.na(clusters)))
      stop("'clusters' contains NA values.", call. = FALSE)
  }
  if (!is.null(blocks)) {
    if (length(blocks) != n)
      stop("'blocks' must have the same length as the number of observations.", call. = FALSE)
    if (any(is.na(blocks)))
      stop("'blocks' contains NA values.", call. = FALSE)
  }
  if (paired && !is.null(blocks))
    stop("Cannot use both 'paired' and 'blocks'.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.validate_clusters_treatment <- function(T, clusters) {
  constant <- tapply(T, clusters, function(x) length(unique(x)) == 1L)
  if (!all(constant)) {
    bad <- names(constant)[!constant]
    stop(sprintf("Treatment is not constant within cluster(s): %s",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
}

#' @keywords internal
#' @noRd
.permute_treatment <- function(T, clusters = NULL, blocks = NULL) {
  n <- length(T)
  original_class <- class(T)
  original_mode <- storage.mode(T)

  if (is.null(clusters) && is.null(blocks)) {
    return(T[sample(n)])
  }
  if (!is.null(clusters) && is.null(blocks)) {
    # Shuffle cluster-level labels, map back
    uclust <- unique(clusters)
    clust_labels <- as.vector(tapply(T, clusters, function(x) x[1L]))
    names(clust_labels) <- uclust
    perm_labels <- clust_labels[sample(length(clust_labels))]
    names(perm_labels) <- uclust
    result <- as.vector(unname(perm_labels[as.character(clusters)]))
    storage.mode(result) <- original_mode
    return(result)
  }
  if (is.null(clusters) && !is.null(blocks)) {
    # Within-block permutation of individual units
    newT <- T
    for (b in unique(blocks)) {
      idx <- which(blocks == b)
      newT[idx] <- T[idx[sample(length(idx))]]
    }
    return(newT)
  }
  # Both clusters and blocks
  newT <- T
  for (b in unique(blocks)) {
    b_idx <- which(blocks == b)
    b_clusters <- unique(clusters[b_idx])
    clust_labels <- as.vector(tapply(T[b_idx], clusters[b_idx], function(x) x[1L]))
    names(clust_labels) <- b_clusters
    perm_labels <- clust_labels[sample(length(clust_labels))]
    names(perm_labels) <- b_clusters
    newT[b_idx] <- as.vector(unname(perm_labels[as.character(clusters[b_idx])]))
  }
  return(newT)
}
