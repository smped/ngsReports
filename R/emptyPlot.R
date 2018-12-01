#' @title Create an empty plot with supplied text
#'
#' @description Create an empty plot with supplied text
#'
#' @details Create plot using \code{theme_void} and only with the supplied text
#'
#' @return A ggplot2 object
#'
#' @examples
#' ngsReports:::.emptyPlot("This is an empty plot")
#'
#' @import ggplot2
#'
#' @keywords internal
#'
.emptyPlot <- function(x){
    ggplot() +
        geom_text(aes(x = 0.5, y = 0.8, label = x)) +
        theme_void() +
        xlim(c(0, 1)) +
        ylim(c(0, 1))
}
