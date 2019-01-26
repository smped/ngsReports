#' @title Plot the Sequence Length Distribution
#'
#' @description Plot the Sequence Length Distribution across one or more FASTQC
#' reports
#'
#' @details
#' This extracts the Sequence Length Distribution from the supplied object and
#' generates a ggplot2 object, with a set of minimal defaults.
#' The output of this function can be further modified using the standard
#' ggplot2 methods.
#' For example, preset axis limits can also be overwritten easily by adding a
#' call to \code{scale_y_continuous} after the call to
#' \code{plotSequenceLengthDistribution}.
#'
#' An alternative interactive plot is available by setting the argument
#' \code{usePlotly = TRUE}.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList},
#' \code{FastqcData}, \code{FastqcDataList} or file path
#' @param usePlotly \code{logical}. Output as ggplot2 or plotly object.
#' @param plotType \code{character}. Can only take the values
#' \code{plotType = "heatmap"} \code{plotType = "line"} or
#' \code{plotType = "cumulative"}
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param counts \code{logical} Should distributions be shown as counts or
#' frequencies (percentages)
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} and
#' \code{usePlotly} are \code{FALSE}. If both \code{cluster} and
#' \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param ... Used to pass additional attributes to theme()
#' @param expand.x Passed to \code{scale_x_discrete}
#' @param heatCol The colour scheme for the heatmap
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Plot as a frequency plot using lines
#' plotSequenceLengthDistribution(fdl)
#'
#' # Or plot the cumulative value
#' plotSequenceLengthDistribution(fdl, plotType = "cumulative")
#'
#' @importFrom dplyr vars
#' @importFrom plotly ggplotly
#' @importFrom plotly layout
#' @importFrom plotly subplot
#' @importFrom viridisLite inferno
#' @import ggplot2
#'
#' @name plotSequenceLengthDistribution
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setGeneric("plotSequenceLengthDistribution", function(
    x, usePlotly = FALSE, labels, ...){
    standardGeneric("plotSequenceLengthDistribution")
}
)
#' @aliases plotSequenceLengthDistribution,character
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "character", function(
    x, usePlotly = FALSE, labels, ...){
    x <- getFastqcData(x)
    plotSequenceLengthDistribution(x, usePlotly, labels, ...)
}
)
#' @aliases plotSequenceLengthDistribution,FastqcFile
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcFile", function(
    x, usePlotly = FALSE, labels, ...){
    x <- getFastqcData(x)
    plotSequenceLengthDistribution(x, usePlotly, labels, ...)
}
)
#' @aliases plotSequenceLengthDistribution,FastqcFileList
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod(
    "plotSequenceLengthDistribution", signature = "FastqcFileList",
    function(x, usePlotly = FALSE, labels, ...){
        x <- getFastqcData(x)
        plotSequenceLengthDistribution(x, usePlotly, labels, ...)
    }
)
#' @aliases plotSequenceLengthDistribution,FastqcData
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod("plotSequenceLengthDistribution", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, plotType = c("line", "cumulative"), ...,
    expand.x = c(0,0.2)){

    df <- Sequence_Length_Distribution(x)
    plotType <- match.arg(plotType)

    if (!length(df)) {
        lenPlot <- .emptyPlot("No Sequence Length Module Detected")
        if (usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
        return(lenPlot)
    }

    labels <- .makeLabels(df, labels, ...)
    df$Filename <- labels[df$Filename]

    ## Add zero counts for lengths either side of the included range
    ## This is only required if a single value exists
    if (nrow(df) == 1) {
        df <- dplyr::bind_rows(
            df,
            dplyr::mutate(df, Lower = Lower - 1, Count = 0),
            dplyr::mutate(df, Lower = Lower + 1, Count = 0)
        )
    }

    df$Lower <- as.integer(df$Lower)
    df <- dplyr::arrange_at(df, vars("Lower"))
    df <- df[c("Filename", "Length", "Lower", "Count")]
    df$Cumulative <- cumsum(df$Count)
    df$Length <- factor(df$Lower, levels = unique(df$Lower))

    ## Sort out some plotting parameters
    xLab <- "Sequence Length (bp)"
    yLab <- c(cumulative = "Cumulative Count", line = "Count")[plotType]
    plotY <- c(cumulative = "Cumulative", line = "Count")[plotType]

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(ggplot2::theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    lenPlot <- ggplot(
        df,
        aes_string("Length", plotY, colour = "Filename", group = "Filename")
    ) +
        geom_line() +
        facet_wrap(~Filename) +
        labs(x = xLab, y = yLab) +
        scale_x_discrete(expand = expand.x) +
        scale_y_continuous(labels = scales::comma) +
        theme_bw() +
        theme(legend.position = "none")

    if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

    if (usePlotly) {
        lenPlot <-
            suppressMessages(plotly::ggplotly(lenPlot, tooltip = c("x", "y")))

        lenPlot <- suppressMessages(
            suppressWarnings(
                plotly::subplot(
                    plotly::plotly_empty(),
                    lenPlot,
                    widths = c(0.14,0.86))
            ))
        lenPlot <- plotly::layout(
            lenPlot,
            xaxis2 = list(title = xLab),
            yaxis2 = list(title = yLab)
        )

    }

    lenPlot

}
)
#' @aliases plotSequenceLengthDistribution,FastqcDataList
#' @rdname plotSequenceLengthDistribution-methods
#' @export
setMethod(
    "plotSequenceLengthDistribution", signature = "FastqcDataList",
    function(
        x, usePlotly = FALSE, labels, counts = FALSE,
        plotType = c("heatmap", "line", "cumulative"), cluster = FALSE,
        dendrogram = FALSE, ..., expand.x = c(0, 0.2), heatCol = inferno(50)){

        df <- Sequence_Length_Distribution(x)

        if (!length(df)) {
            lenPlot <- .emptyPlot("No Sequence Length Module Detected")
            if (usePlotly) lenPlot <- ggplotly(lenPlot, tooltip = "")
            return(lenPlot)
        }

        ## Check for valid plotType
        plotType <- match.arg(plotType)

        ## Drop the suffix, or check the alternate labels
        labels <- .makeLabels(df, labels, ...)

        ## Lengths will probably be binned so define the bins then expand
        ## the range at the lower and upper limits to add zero.
        ## This will enable replication of the default FastQC plot
        ## In reality, this will only be required when there are <2 bins
        lenBins <- stringr::str_sort(unique(df$Length), numeric = TRUE)
        lwr <- upr <- c()
        if (length(lenBins) < 3) {
            lwr <- sprintf("<%i", min(df$Lower))
            upr <- sprintf("%i+", max(df$Upper))
            lwrDf <- tibble(
                Filename = names(labels),
                Length = lwr,
                Count = 0
            )
            uprDf <- tibble(
                Filename = names(labels),
                Length = upr,
                Count = 0
            )
            df <- bind_rows(df, lwrDf, uprDf)
        }

        ## Now spread the gather to fill zeros in any missing bins
        df <- df[c("Filename", "Length", "Count")]
        df <- tidyr::spread(df, "Length", "Count", fill = 0)
        df <- tidyr::gather(df, "Length", "Count", -Filename)

        ## Sort by length bins
        df$Length <- factor(df$Length, levels = c(lwr, lenBins, upr))
        df <- droplevels(df)
        df <- dplyr::arrange(df, Filename, Length)

        ## Get the cumulative count
        df <- dplyr::group_by(df, Filename)
        df <- dplyr::mutate(df, Cumulative = cumsum(Count))
        if (!counts) {
            df <- dplyr::mutate(
                df,
                Cumulative = Cumulative / max(Cumulative),
                Freq = Count / sum(Count))
        }
        df <- dplyr::ungroup(df)
        ## Round the values for better plotting
        df <- dplyr::mutate_if(df, is.double, round, 4)

        ## Get any arguments for dotArgs that have been set manually
        dotArgs <- list(...)
        allowed <- names(formals(ggplot2::theme))
        keepArgs <- which(names(dotArgs) %in% allowed)
        userTheme <- c()
        if (length(keepArgs) > 0)
            userTheme <- do.call(theme, dotArgs[keepArgs])
        ## Rotate labels if >3 lengths
        rot <- ifelse(length(lenBins) > 1, 90, 0)

        if (plotType %in% c("line", "cumulative")) {

            ## Decide whether to plot the Count or cumulative sum
            ## and set all labels
            plotY <- dplyr::case_when(
                plotType == "cumulative" ~ "Cumulative",
                plotType == "line" & counts ~ "Count",
                plotType == "line" & !counts ~ "Freq"
            )
            yLab <- dplyr::case_when(
                plotType == "cumulative" & counts ~ "Cumulative Count",
                plotType == "cumulative" & !counts ~ "Cumulative (%)",
                plotType == "line" & counts ~ "Count",
                plotType == "line" & !counts ~ "Percent (%)"
            )
            yLabelFun <- ifelse(counts, scales::comma, scales::percent)

            df$Filename <- labels[df$Filename]
            lenPlot <- ggplot(
                df,
                aes_string(
                    x = "Length",
                    y = plotY,
                    colour = "Filename",
                    group = "Filename")
            ) +
                geom_line() +
                labs(y = yLab) +
                scale_y_continuous(labels = yLabelFun) +
                theme_bw()
        }

        if (plotType == "heatmap") {

            if (dendrogram && !cluster) {
                message("cluster will be set to TRUE when dendrogram = TRUE")
                cluster <- TRUE
            }
            plotVal <- ifelse(counts, "Count", "Freq")
            fillLab <- ifelse(counts, "Count", "Percent (%)")
            fillLabelFun <- ifelse(counts, scales::comma, scales::percent)

            ## Now define the order for a dendrogram if required
            ## This only applies to a heatmap
            key <- names(labels)
            if (cluster) {
                cols <- c("Filename", "Length", plotVal)
                clusterDend <-
                    .makeDendrogram(df[cols], "Filename","Length", plotVal)
                key <- labels(clusterDend)
            }

            ## Now set everything as factors
            df$Filename <- factor(labels[df$Filename], levels = labels[key])

            ## The additional bins are not really required for a heatmap
            df <- subset(df, subset = Length %in% lenBins)
            df <- droplevels(df)

            lenPlot <- ggplot(
                df,
                aes_string("Length", "Filename", fill = plotVal)
            ) +
                geom_tile() +
                labs(fill = fillLab) +
                scale_fill_gradientn(colours = heatCol, labels = fillLabelFun) +
                scale_y_discrete(labels = labels, expand = c(0, 0))
        }

        ## Set theme elements which are common to all plot types
        lenPlot <- lenPlot +
            scale_x_discrete(expand = expand.x) +
            labs(x = "Sequence Length") +
            theme(
                axis.text.x = element_text(angle = rot, hjust = 1, vjust = 0.5)
            )

        if (!is.null(userTheme)) lenPlot <- lenPlot + userTheme

        if (usePlotly) {

            ## Hide the legend
            lenPlot <- lenPlot + theme(legend.position = "none")

            if (plotType %in% c("line", "cumulative")) {

                ttip <- c("x", "y", "colour")
                lenPlot <- suppressMessages(
                    suppressWarnings(
                        plotly::ggplotly(lenPlot, tooltip = ttip)
                    )
                )
            }
            else{

                pwfCols <- pwf
                lenPlot <- lenPlot  +
                    theme(
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        panel.background = element_blank()
                    )

                status <- getSummary(x)
                status <-
                    subset(status, Category == "Sequence Length Distribution")
                status$Filename <- labels[status$Filename]
                status$Filename <-
                    factor(status$Filename, levels = levels(df$Filename))
                status <- dplyr::right_join(
                    status, unique(df["Filename"]), by = "Filename"
                )

                ## Draw the PWF status as a sideBar
                sideBar <- .makeSidebar(status, key, pwfCols)

                ## Plot dendrogram component
                if (dendrogram) {
                    dx <- ggdendro::dendro_data(clusterDend)
                    dendro <- .renderDendro(dx$segments)
                }
                else{
                    dendro <- plotly::plotly_empty()
                }

                ## Layout the final plot
                lenPlot <- suppressMessages(
                    suppressWarnings(
                        plotly::subplot(
                            dendro, sideBar, lenPlot,
                            widths = c(0.1, 0.08, 0.82),
                            margin = 0.001,
                            shareY = TRUE,
                            titleX = TRUE)
                    )
                )
            }

        }
        lenPlot
    }
)
