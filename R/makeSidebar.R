#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom plotly ggplotly
#'
makeSidebar <- function(status, key, pwfCols){
  nx <- length(status$Filename)
  # Make the basic plot
  sideBar <- ggplot(status, aes(x = 1, y = Filename, key = key)) +
    geom_tile(aes_string(fill = "Status")) +
    geom_hline(yintercept = seq(1.5, nx), colour = "grey20", size = 0.2) +
    scale_fill_manual(values = getColours(pwfCols)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position="none",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
  # Convert to plotly
  suppressWarnings(
    suppressMessages(ggplotly(sideBar, tooltip = c("Status", "Filename")))
    )
}
