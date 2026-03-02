# Internal helpers shared across plot functions

# Suppress R CMD check notes
utils::globalVariables(character(0))

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
