#' @title  A commonly used (hidden) function for setting up dendrograms for interactive plots
#'
#' @description A commonly used (hidden) function for setting up dendrograms for interactive plots
#'
#' @details Create plot using \code{theme_dendro}
#'
#' @param df A `data.frame` as required
#' 
#' @return A ggplot2 object
#'
#' @import ggplot2
#'
#' @keywords internal
#'
ggdend <- function(df) {
    ggplot() +
        geom_segment(data = df, aes_string("x","y", xend = "xend", yend = "yend")) +
        ggdendro::theme_dendro()
}
