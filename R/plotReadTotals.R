#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' The raw data from \code{\link{readTotals}} can otherwise be used to manually
#' create a plot.
#'
#' Duplication levels are based on the value shown on FASTQC reports at the
#' top of the DeDuplicatedTotals plot, which is known to be inaccurate.
#' As it still gives a good guide as to sequence diversity it is included as
#' the default. This can be turned off by setting \code{duplicated = FALSE}.
#'
#' @param x Can be a \code{FastqcData}, \code{FastqcDataList} or file paths
#' @param usePlotly \code{logical} Default \code{FALSE} will render using
#' ggplot. If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param duplicated \code{logical}. Include deduplicated read total estimates
#' to plot charts
#' @param bars If \code{duplicated = TRUE}, show unique and deduplicated reads
#' as "stacked" or "adjacent".
#' @param barCols Colours for duplicated and unique reads.
#' @param expand.x Output from \code{expand_scale()} controlling x-axis
#' expansion. Alternatively can be a numeric vector of length 4
#' @param ... Used to pass additional attributes to theme()
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
#' # Plot the Read Totals showing estimated duplicates
#' plotReadTotals(fdl)
#'
#' # Plot the Read Totals without estimated duplicates
#' plotReadTotals(fdl, duplicated = FALSE)
#'
#' @return Returns a ggplot or plotly object
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#'
#'
#' @name plotReadTotals
#' @rdname plotReadTotals-methods
#' @export
setGeneric("plotReadTotals", function(
    x, usePlotly = FALSE, labels, duplicated = TRUE,
    bars = c("stacked", "adjacent"), barCols = c("red","blue"),
    expand.x = expand_scale(mult = c(0, 0.02)), ...){
    standardGeneric("plotReadTotals")
}
)
#' @aliases plotReadTotals,ANY
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "ANY", function(
    x, usePlotly = FALSE, labels, duplicated = TRUE,
    bars = c("stacked", "adjacent"), barCols = c("red","blue"),
    expand.x = expand_scale(mult = c(0, 0.02)), ...){

    if (length(x) == 1)
        stop("plotReadTotals cannot be called on a single .FastqcFile")

    x <- FastqcDataList(x)
    plotReadTotals(
        x, usePlotly, labels, duplicated, bars,  barCols, expand.x, ...)
}
)
#' @aliases plotReadTotals,FastqcDataList
#' @rdname plotReadTotals-methods
#' @export
setMethod("plotReadTotals", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, duplicated = TRUE,
    bars = c("stacked", "adjacent"), barCols = c("red","blue"),
    expand.x = expand_scale(mult = c(0, 0.02)), ...){

    stopifnot(is.logical(duplicated))
    df <- readTotals(x)

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, ...)
    df$Filename <- labels[df$Filename]
    df$Filename <- factor(df$Filename, levels = unique(df$Filename))

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    ## Get the colours for the barplot
    barCols <- tryCatch(barCols[seq_len(duplicated + 1)])

    ## Check the axis expansion
    stopifnot(length(expand.x) == 4, is.numeric(expand.x))
    xMax <- max(df$Total_Sequences)
    xLab <- "Read Totals"

    if (!duplicated) {

        rtPlot <- ggplot(df, aes_string("Filename", "Total_Sequences")) +
            geom_bar(stat = "identity", fill = barCols) +
            scale_y_continuous(
                labels = scales::comma,
                limits = c(0, xMax),
                expand = expand.x
            ) +
            labs(y = xLab) +
            coord_flip() +
            theme_bw()

        if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme
        if (usePlotly) rtPlot <- plotly::ggplotly(rtPlot)

    }

    if (duplicated) {

        bars <- match.arg(bars)

        ## Add the information to a joined data.frame
        deDup <- getModule(x, "Total_Deduplicated_Percentage")
        deDup$Filename <- labels[deDup$Filename]
        deDup$Filename <- factor(deDup$Filename, levels = levels(df$Filename))
        names(deDup) <- gsub("Total_Deduplicated_", "", names(deDup))
        df <- dplyr::left_join(deDup, df, by = "Filename")

        ##Setup the df for plotting
        types <- c("Unique", "Duplicated")
        df$Unique <- df$Percentage*df$Total_Sequences/100
        df$Unique <- round(df$Unique, 0)
        df$Duplicated <- df$Total_Sequences - df$Unique
        df <- df[c("Filename", types)]
        df <- tidyr::gather(df, "Type", "Total", tidyselect::one_of(types))

        barPos <- c(adjacent = "dodge", stacked = "stack")[bars]

        ## The x-axis expansion needs to be reset for this one
        if (bars == "adjacent") xMax <- max(df$Total)*(1 + expand.x[[1]])

        ## Make the plot
        rtPlot <- ggplot(df, aes_string("Filename", "Total", fill = "Type")) +
            geom_bar(stat = "identity", position = barPos) +
            scale_y_continuous(
                labels = scales::comma,
                limits = c(0, xMax),
                expand = expand.x
            ) +
            scale_fill_manual(values = barCols) +
            labs(y = xLab) +
            coord_flip() +
            theme_bw()

        ## Add common themes & labels
        if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

        if (usePlotly) {

            # Hide the legend
            rtPlot <- rtPlot + theme(legend.position = "none")
            # Render as a plotly object
            rtPlot <- suppressMessages(
                suppressWarnings(plotly::ggplotly(rtPlot))
            )
        }

    }

    ## Draw the plot
    rtPlot
}
)
