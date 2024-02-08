#' @title Plot the PASS/WARN/FAIL information
#'
#' @description Extract the PASS/WARN/FAIL summaries and plot them
#'
#' @details This uses the standard ggplot2 syntax to create a three colour plot.
#' The output of this function can be further modified using the standard
#' ggplot2 methods if required.
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or character
#' vector of file paths
#' @param pwfCols Object of class [PwfCols()] containing the colours
#' for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly `logical`. Generate an interactive plot using plotly
#' @param cluster `logical` default `FALSE`. If set to `TRUE`,
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram `logical` redundant if `cluster` is `FALSE`
#' if both `cluster` and `dendrogram` are specified as `TRUE`
#' then the dendrogram will be displayed.
#' @param ... Used to pass various potting parameters to theme.
#' @param gridlineWidth,gridlineCol Passed to geom_hline and geom_vline to
#' determine width and colour of gridlines
#'
#' @return A ggplot2 object (`usePlotly = FALSE`)
#' or an interactive plotly object (`usePlotly = TRUE`)
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Check the overall PASS/WARN/FAIL status
#' plotSummary(fdl)
#'
#' @docType methods
#'
#' @importFrom dplyr bind_cols
#' @importFrom grid unit
#'
#' @name plotSummary
#' @rdname plotSummary-methods
#' @export
setGeneric("plotSummary", function(
        x, usePlotly = FALSE, labels, pwfCols, cluster = FALSE, dendrogram = FALSE,
        ...){
    standardGeneric("plotSummary")
}
)
#' @rdname plotSummary-methods
#' @export
setMethod("plotSummary", signature = "ANY", function(
        x, usePlotly = FALSE, labels, pwfCols, cluster = FALSE, dendrogram = FALSE,
        ...){
    .errNotImp(x)
}
)
#' @rdname plotSummary-methods
#' @export
setMethod("plotSummary", signature = "character", function(
        x, usePlotly = FALSE, labels, pwfCols, cluster = FALSE, dendrogram = FALSE,
        ...){
    if (length(x) == 1)
        stop("plotSummary() can only be called on two or more files.")
    x <- FastqcDataList(x)
    plotSummary(x, usePlotly, labels, pwfCols, cluster, dendrogram, ...)
}
)
#' @rdname plotSummary-methods
#' @export
setMethod("plotSummary", signature = "FastqcDataList", function(
        x, usePlotly = FALSE, labels, pwfCols, cluster = FALSE,
        dendrogram = FALSE, ..., gridlineWidth = 0.2, gridlineCol = "grey20"
){

    df <- getSummary(x)

    if (missing(pwfCols)) pwfCols <- ngsReports::pwf
    stopifnot(.isValidPwf(pwfCols))
    fillCol <- getColours(pwfCols)

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(x, labels, ...)
    labels <- labels[names(labels) %in% df$Filename]

    ## Set factor levels
    df$Category <- factor(df$Category, levels = unique(df$Category))
    df$Status <- factor(df$Status, levels = rev(c("PASS", "WARN", "FAIL")))
    df$StatusNum <- as.integer(df$Status)

    ## Make sure cluster is TRUE if the dendrogram is requested
    if (dendrogram && !cluster) {
        message("cluster will be set to TRUE when dendrogram = TRUE")
        cluster <- TRUE
    }

    ## Now define the order for a dendrogram if required
    ## This only applies to a heatmap
    key <- names(labels)
    if (cluster) {
        cols <- c("Filename", "Category", "StatusNum")
        clusterDend <- .makeDendro(
            df[cols], "Filename", "Category", "StatusNum"
        )
        key <- labels(clusterDend)
    }

    ## Set factor levels accordingly
    df$Filename <- factor(labels[df$Filename], levels = labels[key])

    ## Create the basic plot
    sumPlot <- ggplot(df, aes(Category, Filename, fill = Status)) +
        geom_tile(colour = gridlineCol, linewidth = gridlineWidth) +
        scale_fill_manual(values = fillCol) +
        labs(x = "QC Category", y = "Filename") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw() +
        theme(axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5
        ))
    sumPlot <- .updateThemeFromDots(sumPlot, ...)

    if (usePlotly) {

        ## Set the dimensions of the plot
        ny <- length(x)
        nx <- length(unique(df$Category))
        mar <- unit(c(0.01, 0.01, 0.01, 0.04), "npc")
        cmn_thm <- theme(
            axis.title = element_blank(), plot.margin = mar,
            legend.position = "none"
        )

        if (dendrogram) {

            ## Remove the legend and labels for plotly
            sumPlot <- sumPlot +
                cmn_thm +
                theme(
                    axis.text.y = element_blank(), axis.ticks = element_blank()
                )

            ## Get the dendrogram sorted out
            dx <- ggdendro::dendro_data(clusterDend)
            dendro <- .renderDendro(dx$segments)

            ## Now layout the plotly version
            sumPlot <- suppressWarnings(
                suppressMessages(
                    plotly::subplot(
                        dendro, sumPlot, widths = c(0.1, 0.9), margin = 0.001,
                        shareY = TRUE
                    )
                )
            )
        }
        else{

            ## Remove the legend and set the margin for plotly
            ## Generate as a single plot
            sumPlot <- plotly::ggplotly(sumPlot + cmn_thm)

        }

    }

    sumPlot

}
)
