#' @title Shortcut for making the status sidebar
#'
#' @description Shortcut for making the status sidebar
#'
#' @param status A data.frame with columns 'Filename' & 'Status'
#' @param key A vector of values corresponding to the Filename column
#' @param pwfCols An object of class PwfCols
#'
#' @return
#' A plotly object. The sidebar for an interactive plot showing PASS/WARN/FAIL
#' status for each file.
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#'
#' @keywords internal
#'
.makeSidebar <- function(status, key, pwfCols){

    stopifnot(.isValidPwf(pwfCols))
    nx <- length(status$Filename)
    ## make sure status is in right order so key can apply
    ## This only works because the factor levels of the 'Filename' column
    ## correspond to the order of the key as determined earlier the plotting
    ## functions. This step is essentially redundant
    status <- status[order(status$Filename),]
    ## Make the basic plot
    sideBar <-
        ggplot(status, aes_string(x = "1", y = "Filename", key = "key")) +
        geom_tile(aes_string(fill = "Status")) +
        geom_hline(yintercept = seq(1.5, nx), colour = "grey20", size = 0.2) +
        scale_fill_manual(values = getColours(pwfCols)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme(
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()
        )

    ## Convert to plotly
    suppressWarnings(
        suppressMessages(ggplotly(sideBar, tooltip = c("y", "fill")))
    )
}
