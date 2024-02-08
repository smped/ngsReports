#' @title Split elements of a vector into a data.frame
#'
#' @description Split elements of a character vector by the tab separator
#'
#' @details This will split a vector into a data.frame checking that every line
#' has the same number of separators.
#' By default the first element will be set as the column names.
#'
#' This is designed to take input from `readLines()`
#'
#' @param x A character vector
#' @param firstRowToNames logical Should the first element be used for column
#' names
#' @param tab character The string used to represent the tab symbol
#'
#' @return
#' A data frame
#'
#' @examples
#' x <- c("ColA\tColB", "Value1\tValue2")
#' ngsReports:::.splitByTab(x, firstRowToNames = TRUE)
#' ngsReports:::.splitByTab(x, firstRowToNames = FALSE)
#'
#' @keywords internal
#'
.splitByTab <- function(x, firstRowToNames = TRUE, tab = "\\t"){

    stopifnot(is.character(x))

    ## Check for the tab marker in every line
    linesWithTab <- stringr::str_detect(x, tab)
    if (sum(linesWithTab) != length(x))
        stop("Some elements of x are missing the tab separator")

    ## Take the first element as defining the number of columns
    nCol <- stringr::str_count(x[1], pattern = tab) + 1

    ## Count the number of tabs in each line
    nTabs <- stringr::str_count(string = x, pattern = tab)
    if (any(nTabs != (nCol - 1)))
        stop("Differing number of delimiters in some rows")

    if (firstRowToNames) {

        ## Get the first element as a vector of names
        nm <- stringr::str_split_fixed(x[1], pattern = tab, n = nCol)

        ## Split the remainder
        df <- stringr::str_split_fixed(x[-1], pattern = tab, n = nCol)
        colnames(df) <- nm
    }
    else {
        df <- stringr::str_split_fixed(x, pattern = tab, n = nCol)
    }
    ## Return a generic data.frame
    ## This leaves tidying to each module
    as.data.frame(df, stringsAsFactors = FALSE)
}

#' @title Add a percentage sign
#'
#' @description Add a percentage sign to the end of a string
#'
#' @param x Any vector
#'
#' @return character vector
#'
#' @examples
#'
#' x <- 1:10
#' ngsReports:::.addPercent(x)
#'
#' @keywords internal
#'
.addPercent <- function(x){
    if (is.factor(x)) message("Factors will be converted to characters")
    paste0(x, "%")
}

#' @title Create an empty plot with supplied text
#'
#' @description Create an empty plot with supplied text
#'
#' @details Create plot using `theme_void` and only with the supplied text
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

#' @title Make the dendrogram for heatmap-style plots
#'
#' @description Set the clusters for heatmap-style interactive plots
#'
#' @param df The data frame to be clustered
#' @param rowVal The rows to be clustered
#' @param colVal The value which will become column names
#' @param value The value to use for the clustering
#'
#' @return A dendrogram
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#' cols <- c("Filename", "Position", "Illumina_Universal_Adapter")
#' ac <- getModule(fileList, "Adapter_Content")[cols]
#' ngsReports:::.makeDendro(df = ac,
#'                             rowVal = "Filename",
#'                             colVal = "Position",
#'                             value = "Illumina_Universal_Adapter")
#'
#' @importFrom stats as.formula
#'
#' @keywords internal
.makeDendro <- function(df, rowVal, colVal, value){

    cols <- c(rowVal, colVal, value)
    stopifnot(all(cols %in% names(df)))
    df <- df[cols]
    fm <- as.formula(paste0("`", rowVal, "`~`", colVal, "`"))
    mat <- reshape2::acast(df, fm, value.var = value)
    mat[is.na(mat)] <- 0
    clust <- hclust(dist(mat), method = "ward.D2")
    as.dendrogram(clust)

}

#' @title Perform the checks and return the labels
#'
#' @description Checks for the presence of labels and returns defaults
#'
#' @details Takes a named vector of labels and checks for the correct fields.
#' If no vector is supplied, returns the file names missing the specified
#' pattern, which defaults to removing the suffixes fastq(.gz), fq(.gz),
#' bam, sam or cram.
#'
#' @param x A data.frame with a column titled "Filename"
#' @param labels Named vector of labels for plotting
#' @param pattern character Regular expression to remove from filenames
#' @param col character Column to use for generating labels
#' @param ... Not used
#'
#' @return Named character vector
#'
#' @examples
#' f <- paste0(c("File1", "File2"), ".fastq")
#' df <- data.frame(Filename = f, stringsAsFactors = FALSE)
#' ngsReports:::.makeLabels(df)
#'
#' @keywords internal
.makeLabels <- function(
        x, labels, pattern = ".(fast|fq|bam|sam|cram).*", col = "Filename", ...
){

    if (is(x, "FastqcDataList") | is(x, "FastqcData")) {
        ## Form a single column data.frame
        x <- structure(
            list(fqName(x)), names = col, row.names = seq_along(x),
            class = "data.frame"
        )
    }

    stopifnot(is(x, "data.frame"))

    col <- match.arg(col, colnames(x))

    ## If no labels are provided, just remove the file suffix as
    ## determined by the supplied pattern
    if (missing(labels)) {
        labels <- structure(
            gsub(pattern, "", unique(x[[col]])), # Remove the pattern
            names = unique(x[[col]]) # Ensure a named vector
        )
    }
    if (!all(x[[col]] %in% names(labels)))
        stop("Names of supplied labels must match all filenames.")
    if (any(duplicated(labels))) stop("Labels must be unique.")

    ## Now return only the supplied labels which are in the df
    labels[names(labels) %in% x[[col]]]
}

#' @title Shortcut for making the status sidebar
#'
#' @description Shortcut for making the status sidebar
#'
#' @param status A data.frame with columns 'Filename' & 'Status'
#' @param key A vector of values corresponding to the Filename column
#' @param pwfCols An object of class PwfCols
#' @param usePlotly If TRUE, output is a plotly panel
#'
#' @return
#' if `usePlotly = TRUE`, a plotly object. The sidebar for an interactive plot
#' showing PASS/WARN/FAIL status for each file. If `usePlotly = FALSE` the
#' underlying `ggplot` object will be returned.
#'
#' @import ggplot2
#'
#' @keywords internal
#'
.makeSidebar <- function(status, key, pwfCols, usePlotly = TRUE){

    stopifnot(.isValidPwf(pwfCols))
    nx <- length(status$Filename)
    ## make sure status is in right order so key can apply
    ## This only works because the factor levels of the 'Filename' column
    ## correspond to the order of the key as determined earlier the plotting
    ## functions. This step is now essentially redundant
    status <- status[order(status$Filename),]
    ## Make the basic plot
    sideBar <- ggplot(status, aes(1, Filename, key = key)) +
        geom_tile(aes(fill = Status)) +
        geom_hline(yintercept = seq(1.5, nx), colour = "grey20", linewidth = 0.2) +
        scale_fill_manual(values = getColours(pwfCols)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        theme(
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none", axis.title = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank()
        )

    ## Convert to plotly
    if (usePlotly) {
        sideBar <- suppressWarnings(
            suppressMessages(plotly::ggplotly(sideBar, tooltip = c("y", "fill")))
        )
    }
    sideBar
}

#' @title  Set up dendrograms for interactive plots
#'
#' @description A commonly used (hidden) function for setting up dendrograms
#' for interactive plots. based on code found at
#' https://plot.ly/ggplot2/ggdendro-dendrograms/
#'
#' @details Create plot using `theme_dendro`
#'
#' @param df A `data.frame` as required
#'
#' @return A plotly object
#'
#' @import ggplot2
#' @importFrom ggdendro theme_dendro
#'
#' @keywords internal
#'
.renderDendro <- function(df) {
    ## Based on the example ggdend
    dendro <- ggplot() +
        geom_segment(data = df, aes(x, y, xend = xend, yend = yend)) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0.5)) +
        theme_dendro()
    plotly::ggplotly(dendro, tooltip = NULL)
}


#' @title Construct a gradient using PwfCols
#'
#' @description Construct a gradient using PwfCols
#'
#' @details This constructs a list of arguments for passing to
#' `scale_fill_gradientn()` using the values contained in the
#' fill aesthetic and the supplied breaks for PASS/WARN/FAIL criteria.
#'
#' @param vals The values which need to have the scale generated
#' @param pwfCols An object of class PwfCols
#' @param breaks The breaks for the PWF bins
#' @param passLow Is the PASS category at the low or high end of the numeric
#' range
#' @param na.value The colour to plot for missing values
#'
#' @return
#' Returns a ggplot list
#'
#' @include PwfCols.R
#' @importFrom grDevices colorRampPalette
#'
#' @keywords internal
#'
.makePwfGradient <- function(
        vals, pwfCols, breaks = c(0, 5, 10, 100),  passLow = TRUE,
        na.value = "white"
){

    ## passLow defines whether pass is the low score or the high score
    ## organise the colours based on this
    o <- seq_len(4)
    if (!passLow) o <- rev(o)
    gradCols <- getColours(pwfCols)[o] # Get the default gradient colours

    ## Find which of the pwf bins are present
    bins <- cut(vals, breaks = breaks, include.lowest = TRUE)
    bins <- range(as.integer(bins))
    bins <- unique(bins)

    ## Create an even sequence between the min & max of the range
    n <- seq(breaks[min(bins)], breaks[max(bins) + 1], length.out = 101)
    n <- cut(n, breaks = breaks)
    n <- split(n, n)

    ## Now create a colour vector between the extreme points
    cols <- lapply(bins, function(x){
        l <- length(n[[x]]) + 1
        colorRampPalette(gradCols[c(x,x + 1)])(l)[-l]
    })
    cols <- as.character(c(unlist(cols), gradCols[max(bins) + 1]))

    ## Remove any breaks outside of the range
    breaks <- breaks[seq(min(bins), max(bins) + 1)]

    ## Return a list for passing to do.call("scale_fill_gradientn", args)
    list(
        colours = cols, breaks = breaks,
        limits = range(breaks),
        na.value = na.value
    )

}

#' @title Prepare the final heatmap for plotting
#'
#' @description Add dendrogarm & status bar to ggplot2 heatmap
#'
#' @param x a ggplot2 heatmap produced by ngsReports
#' @param status a tibble with the columns Filename and Status
#' @param segments a dendrogram produced during clustering of samples
#' @param usePlotly logical(1)
#' @param hv character vector of fields to include in hoverinfo
#'
#' @return
#' Either a ggplot2 object assembled using patchwork, or an interactive plotly
#' object
#'
#' @import patchwork
#' @import ggplot2
#'
#' @keywords internal
.prepHeatmap <- function(
        x, status, segments, usePlotly, heat_w = 8, pwf, hv = NULL
) {

    stopifnot(is(x, "gg"))
    hasStatus <- as.logical(nrow(status))
    stopifnot(all(c("x", "y", "xend", "yend") %in% colnames(segments)))

    if (missing(pwf)) pwf <- ngsReports::pwf

    ## Create the dendrogram if required. This is independent of plotly
    add_dend <- nrow(segments) > 0
    panel_w <- c(1, heat_w)
    if (add_dend) {
        n <- max(segments$xend)
        panel_w <- c(1, 1, heat_w)
        dendPlot <- ggplot(segments) +
            geom_segment(aes(x = yend, y = xend, xend = y, yend = x)) +
            scale_x_reverse(expand = expansion(0)) +
            scale_y_continuous(limits = c(0, n) + 0.5, expand = expansion(0)) +
            labs(x = "", y = c()) +
            theme_minimal() +
            theme(
                panel.grid = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(),
                plot.margin = unit(c(5.5, 0, 5.5, 5.5), "points")
            )
    }

    x_lab <- x$labels$x
    if (hasStatus) {
        stopifnot(all(c("Filename", "Status") %in% colnames(status)))
        ## Now create the sideBar
        sideBar <- .makeSidebar(status, levels(status$Filename), pwf, usePlotly)

        if (!usePlotly) {

            sideBar <- sideBar +
                theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "points"))
            out <- sideBar
            if (add_dend) out <- dendPlot + sideBar
            out <- out + x +  plot_layout(widths = panel_w)

        } else {

            ## Setup additional formatting
            x <- x + theme(
                plot.title = element_text(hjust = 0.5),
                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                legend.position = "none"
            )
            title_x <- 1 - 0.5 * (heat_w + 1) / sum(panel_w)
            panel_w <- panel_w / sum(panel_w)
            if (is.null(hv)) hv <- "all"
            x <- plotly::ggplotly(x, tooltip = hv)

            if (add_dend) {
                out <- suppressWarnings(
                    suppressMessages(
                        plotly::subplot(
                            .renderDendro(segments), sideBar, x,
                            widths = panel_w, margin = 0.001, shareY = TRUE
                        )
                    )
                )
                out <- plotly::layout(
                    out, title = list(x = title_x), xaxis3 = list(title = x_lab),
                    margin = list(b = 50, t = 50)
                )
            } else {
                out <- suppressWarnings(
                    suppressMessages(
                        plotly::subplot(
                            sideBar, x,
                            widths = panel_w, margin = 0.001, shareY = TRUE
                        )
                    )
                )
                out <- plotly::layout(
                    out, title = list(x = title_x),
                    xaxis2 = list(title = x_lab), margin = list(b = 50, t = 50)
                )
            }
        }

    } else {
        panel_w <- panel_w[c(1, 3)]
        if (!usePlotly) {

            out <- x
            if (add_dend) out <- dendPlot + x + plot_layout(widths = panel_w)

        } else {

            ## Setup additional formatting
            x <- x + theme(
                plot.title = element_text(hjust = 0.5),
                axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                legend.position = "none"
            )
            title_x <- 1 - 0.5 * (heat_w + 1) / sum(panel_w)
            panel_w <- panel_w / sum(panel_w)
            if (is.null(hv)) hv <- "all"
            x <- plotly::ggplotly(x, tooltip = hv)

            if (add_dend) {
                out <- suppressWarnings(
                    suppressMessages(
                        plotly::subplot(
                            .renderDendro(segments), x,
                            widths = panel_w, margin = 0.001, shareY = TRUE
                        )
                    )
                )
                out <- plotly::layout(
                    out, title = list(x = title_x),
                    xaxis2 = list(title = x_lab), margin = list(b = 50, t = 50)
                )
            } else {
                out <- x
            }

        }
    }

    out

}

#' Add custom theme elements from dotArgs
#'
#' @param p ggplot object
#' @param ... Standard dot arguments
#'
#' @return ggplot2 object
#' @keywords internal
.updateThemeFromDots <- function(p, ...){
    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
    ## Add the basic customisations
    if (!is.null(userTheme)) p <- p + userTheme
    p
}


#' @title Hide PWF tooltips from line plots
#' @description Hide tooltips from PWF rectangles in line plots
#' @param x plotlyObject$x$data
#' @return plotlyObject$x$data
#' @keywords internal
.hidePWFRects <- function(x){
    ## If there is a name component & it contains
    ## PASS/WARN/FAIL set the hoverinfo to none
    if (any(grepl("(PASS|WARN|FAIL|Status)", x$text))) {
        x$hoverinfo <- "none"
        x$text <- ""
    }
    if (!is.null(x$name)) {
        if (any(grepl("(PASS|WARN|FAIL|Status)", x$name))) {
            x$hoverinfo <- "none"
            x$text <- ""
        }
    }
    x
}
