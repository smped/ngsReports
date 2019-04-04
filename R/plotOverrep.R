#' @title Plot a summary of Over-represented Sequences
#'
#' @description Plot a summary of Over-represented Sequences for a set of
#' FASTQC reports
#'
#' @details Percentages are obtained by simply summing those within a report.
#' Any possible double counting by FastQC is ignored for the purposes of a
#' simple approximation.
#'
#' Plots generated from a \code{FastqcData} object will show the top \code{n}
#' sequences grouped by their predicted source & coloured by whether the
#' individual sequence would cause a WARN/FAIL.
#'
#' Plots generated from a \code{FastqcDataList} group sequences by predicted
#' source and summarise as a percentage of the total reads.
#'
#' @param x Can be a \code{FastqcData}, \code{FastqcDataList} or file paths
#' @param usePlotly \code{logical} Default \code{FALSE} will render using
#' ggplot. If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param n The number of sequences to plot from an individual file
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours
#' for PASS/WARN/FAIL
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE}
#' then the dendrogram will be displayed.
#' @param ... Used to pass additional attributes to theme() and between methods
#' @param expand.x,expand.y Output from \code{expand_scale()} or numeric
#' vectors of length 4. Passed to \code{scale_*_continuous()}
#' @param paletteName Name of the palette for colouring the possible sources
#' of the overrepresented sequences. Must be a palette name from
#' \code{RColorBrewer}
#'
#' @return A standard ggplot2 object
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
#' # Another example which isn't ideal
#' plotOverrep(fdl)
#'
#' @docType methods
#'
#' @importFrom tidyr spread
#' @importFrom plotly layout
#' @importFrom grDevices rgb
#' @import ggplot2
#'
#' @name plotOverrep
#' @rdname plotOverrep-methods
#' @export
setGeneric("plotOverrep", function(
    x, usePlotly = FALSE, labels, pwfCols, ...){
    standardGeneric("plotOverrep")
}
)
#' @rdname plotOverrep-methods
#' @export
setMethod("plotOverrep", signature = "ANY", function(
    x, usePlotly = FALSE, labels, pwfCols, ...){
    .errNotImp(x)
}
)
#' @rdname plotOverrep-methods
#' @export
setMethod("plotOverrep", signature = "character", function(
    x, usePlotly = FALSE, labels, pwfCols, ...){
    x <- FastqcDataList(x)
    if (length(x) == 1) x <- x[[1]]
    plotOverrep(x, usePlotly, labels, pwfCols, ...)
}
)
#' @rdname plotOverrep-methods
#' @export
setMethod("plotOverrep", signature = "FastqcData", function(
    x, usePlotly = FALSE, labels, pwfCols,  n = 10, ...,
    expand.x = expand_scale(mult = c(0, 0.05)),
    expand.y = expand_scale(0, 0.6)){

    df <- getModule(x, "Overrepresented_sequences")

    if (!length(df)) {
        overPlot <- .emptyPlot("No Overrepresented Sequences")
        if (usePlotly) overPlot <- ggplotly(overPlot, tooltip = "")
        return(overPlot)
    }

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, ...)
    df$Filename <- labels[df$Filename]

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    df <- dplyr::top_n(df, n, Percentage)
    df$Status <- cut(
        df$Percentage,
        breaks = c(0, 0.1, 1, 100),
        labels = c("PASS", "WARN", "FAIL")
    )
    df$Possible_Source <-
        gsub(" \\([0-9]*\\% over [0-9]*bp\\)",  "", df$Possible_Source)
    df$Sequence <- factor(df$Sequence, levels = rev(df$Sequence))
    df$Percentage <- round(df$Percentage, 2)
    df <- droplevels(df)

    ## Check the axis expansion
    stopifnot(is.numeric(expand.x), length(expand.x) == 4)
    stopifnot(is.numeric(expand.y), length(expand.y) == 4)

    ## Set plotting parameters
    xLab <- "Percent of Total Reads (%)"
    yLab <- "Overrepresented Sequence"

    ## Sort out the colours & pass/warn/fail breaks
    if (missing(pwfCols)) pwfCols <- getColours(ngsReports::pwf)
    pwfCols <- pwfCols[names(pwfCols) %in% levels(df$Status)]

    overPlot <- ggplot(
        df,
        aes_string(
            x = "Sequence", y = "Percentage", fill = "Status",
            label = "Possible_Source"
        )
    ) +
        geom_bar(stat = "identity") +
        labs(y = xLab, x = yLab) +
        scale_y_continuous(expand = expand.x, labels = .addPercent) +
        scale_x_discrete(expand = expand.y) +
        theme_bw() +
        coord_flip() +
        scale_fill_manual(values = pwfCols)

    ## Only facet is using ggplot. They look bad under plotly
    if (!usePlotly) overPlot <- overPlot +
        facet_grid(Possible_Source~., scales = "free_y", space = "free")

    ## Add the basic customisations
    if (!is.null(userTheme)) overPlot <- overPlot + userTheme

    if (usePlotly) {

        ## Add the customisations for plotly
        overPlot <- overPlot +
            ggtitle(df$Filename[1]) +
            theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.position = "none"
            )

        ## Add the empty plot to align in the shiny app
        overPlot <- suppressMessages(
            suppressWarnings(
                plotly::subplot(
                    plotly::plotly_empty(),
                    overPlot,
                    widths = c(0.14,0.86)
                )
            )
        )
    }

    overPlot

}
)
#' @rdname plotOverrep-methods
#' @export
setMethod("plotOverrep", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, pwfCols, cluster = TRUE, dendrogram = TRUE,
    ..., paletteName = "Set1", expand.x = expand_scale(mult = c(0, 0.05)),
    expand.y = expand_scale(0, 0)){

    df <- getModule(x, "Overrepresented_sequences")

    if (!length(df)) {
        overPlot <- .emptyPlot("No Overrepresented Sequences")
        if (usePlotly) overPlot <- ggplotly(overPlot, tooltip = "")
        return(overPlot)
    }

    if (missing(pwfCols)) pwfCols <- ngsReports::pwf

    ## Drop the suffix, or check the alternate labels
    labels <- .makeLabels(df, labels, ...)

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    Possible_Source <- c() # Here to avoid a NOTE in R CMD check...
    df$Possible_Source <-
        gsub(" \\([0-9]*\\% over [0-9]*bp\\)", "", df$Possible_Source)
    df <- dplyr::group_by(df, Filename, Possible_Source)
    df <- dplyr::summarise(df, Percentage = sum(Percentage))
    df <- dplyr::ungroup(df)
    df$Percentage <- round(df$Percentage, 2)
    lev <- unique(dplyr::arrange(df, Percentage)$Possible_Source)
    df$Possible_Source <- factor(df$Possible_Source, levels = lev)

    if (dendrogram && !cluster) {
        message("cluster will be set to TRUE when dendrogram = TRUE")
        cluster <- TRUE
    }

    ## Now define the order for a dendrogram if required
    key <- names(labels)
    if (cluster) {
        cols <- c("Filename", "Possible_Source", "Percentage")
        clusterDend <-  .makeDendro(
            df[cols], "Filename", "Possible_Source", "Percentage")
        key <- labels(clusterDend)
    }
    ## Now set everything as factors
    df$Filename <- factor(labels[df$Filename], levels = labels[key])
    maxChar <- max(nchar(levels(df$Filename)))

    ## Check the axis expansion
    stopifnot(is.numeric(expand.x), length(expand.x) == 4)
    stopifnot(is.numeric(expand.y), length(expand.y) == 4)

    ## Define the palette
    paletteName <-
        match.arg(paletteName, rownames(RColorBrewer::brewer.pal.info))
    nMax <- RColorBrewer::brewer.pal.info[paletteName, "maxcolors"]
    nSource <- length(levels(df$Possible_Source))
    pal <- RColorBrewer::brewer.pal(nMax, paletteName)
    if (nSource > nMax) {
        pal <- colorRampPalette(pal)(nSource)
    }
    else {
        pal <- pal[seq_len(nSource)]
    }
    names(pal) <- levels(df$Possible_Source)

    ## Set the axis label
    xLab <- "Overrepresented Sequences (% of Total)"

    overPlot <- ggplot(
        df,
        aes_string("Filename", "Percentage", fill = "Possible_Source")
    ) +
        geom_bar(stat = "identity") +
        labs(y = xLab, fill = "Possible Source") +
        scale_y_continuous(expand = expand.x, labels = .addPercent) +
        scale_x_discrete(expand = expand.y) +
        scale_fill_manual(values = pal) +
        theme_bw() +
        coord_flip()

    ## Add the basic customisations
    if (!is.null(userTheme)) overPlot <- overPlot + userTheme

    if (usePlotly) {

        # Remove annotations before sending to plotly
        overPlot <- overPlot +
            theme(
                legend.position = "none",
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank()
            )
        # Prepare the sidebar
        status <- getSummary(x)
        status <- subset(status, Category == "Overrepresented sequences")
        status$Filename <- labels[status$Filename]
        status$Filename <-
            factor(status$Filename, levels = levels(df$Filename))
        status <- dplyr::right_join(
            status,
            dplyr::distinct(df, Filename),
            by = "Filename"
        )
        sideBar <- .makeSidebar(status, key, pwfCols)

        # Prepare the dendrogram
        dendro <- plotly::plotly_empty()
        if (dendrogram) {
            dx <- ggdendro::dendro_data(clusterDend)
            dendro <- .renderDendro(dx$segments)
        }

        # The final interactive plot
        overPlot <- suppressWarnings(
            suppressMessages(
                plotly::subplot(
                    dendro,
                    sideBar,
                    overPlot,
                    margin = 0.001,
                    widths = c(0.08,0.08,0.84)
                )
            )
        )
    }
    overPlot
}
)
