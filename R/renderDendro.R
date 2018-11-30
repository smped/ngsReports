#' @title  Set up dendrograms for interactive plots
#'
#' @description A commonly used (hidden) function for setting up dendrograms for 
#' interactive plots
#'
#' @details Create plot using \code{theme_dendro}
#'
#' @param df A `data.frame` as required
#' 
#' @return A plotly object
#'
#' @import ggplot2
#' @importFrom ggdendro theme_dendro
#' @importFrom plotly ggplotly
#'
#' @keywords internal
#'
renderDendro <- function(df) {
    # Based on the example ggdend
    dendro <- ggplot() +
        geom_segment(data = df, aes_string("x","y", 
                                           xend = "xend", yend = "yend")) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0.5)) +
        theme_dendro() 
    ggplotly(dendro, tooltip = NULL)
}
