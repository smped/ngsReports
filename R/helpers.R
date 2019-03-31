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
#' @param tab character The string used torepresent the tab symbol
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

    stopifnot(setequal(c(rowVal, colVal, value), names(df)))
    fm <- as.formula(paste0(rowVal, "~", colVal))
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
#' @param df A data.frame with a column titled "Filename"
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
    df, labels, pattern = ".(fastq|fq|bam|sam|cram).*", col ="Filename", ...){
    stopifnot(is.data.frame(df))
    col <- match.arg(col, colnames(df))
    ## If no labels are provided, just remove the file suffix as determined by
    ## the supplied pattern
    if (missing(labels)) {
        labels <- structure(
            gsub(pattern, "", unique(df[[col]])), # Remove the pattern
            names = unique(df[[col]]) # Ensure a named vector
        )
    }
    if (!all(df[[col]] %in% names(labels)))
        stop("Names of supplied labels must match all filenames.")
    if (any(duplicated(labels))) stop("Labels must be unique.")
    labels
}

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
    ## functions. This step is now essentially redundant
    status <- status[order(status$Filename),]
    ## Make the basic plot
    sideBar <- ggplot(status, aes_string("1", "Filename", key = "key")) +
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

#' @title  Set up dendrograms for interactive plots
#'
#' @description A commonly used (hidden) function for setting up dendrograms
#' for interactive plots. based on code found at
#' https://plot.ly/ggplot2/ggdendro-dendrograms/
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
.renderDendro <- function(df) {
    ## Based on the example ggdend
    dendro <- ggplot() +
        geom_segment(
            data = df,
            aes_string("x","y", xend = "xend", yend = "yend")
        ) +
        coord_flip() +
        scale_y_reverse(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0.5)) +
        theme_dendro()
    ggplotly(dendro, tooltip = NULL)
}

#' @title Construct a scale using PwfCols
#'
#' @description Construct a scale using PwfCols
#'
#' @details This constructs a ggplot scale using the values contained in the
#' fill aesthetic and the supplied breaks for PASS/WARN/FAIL criteria.
#' As this doesn't follow the conventional ggplot syntax, using more of
#' functional approach, it will be a hidden function.
#'
#' @param vals The values which need to have the scale generated
#' @param pwfCols An object of class PwfCols
#' @param breaks The breaks for the PWF bins
#' @param passLow Is the PASS category at the low or high end of the numeric
#' range
#' @param na.value The colour to plot for missing values
#'
#' @return
#' Returns a ggplot scale object
#'
#' @include PwfCols.R
#' @import ggplot2
#'
#' @keywords internal
#'
.scale_fill_pwf <- function(
    vals, pwfCols, breaks = c(0, 5, 10, 100),  passLow = TRUE,
    na.value = "white"){

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

    ## Return the scale object
    scale_fill_gradientn(
        colours = cols,
        breaks = breaks,
        limits = range(breaks),
        na.value = na.value
    )

}
