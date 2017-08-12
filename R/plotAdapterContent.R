#' @title Draw an Adapter Content Plot
#'
#' @description Draw an Adapter Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the Adapter_Content from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#'
#' Preset axis limits can also be overwritten easily by adding a call to \code{scale_y_continuous}
#' after the call to \code{plotAdapterContent}.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param adapterType A regular expression used to filter which adapter(s) are plotted
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param nc Number of columns to use when faceting by Adapter Type.
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param trimNames \code{logical}. Remove the file suffix from the names displyed in the legend.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param ... Used to pass additional arguments, such as \code{invert} or \code{ignore.case},
#' to \code{grep} when selecting Adapter \code{adapterType}
#'
#' @return A standard ggplot2 object
#'
#' @examples
#'
#' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot
#' plotAdapterContent(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fileName(fdl))
#' plotAdapterContent(fdl, subset = r1)
#'
#' # Plot just the Universal Adapter, and change the y-axis
#' library(ggplot2)
#' plotAdapterContent(fdl, adapterType ="Universal", plotType = "line") +
#' scale_y_continuous()
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#'
#' @export
plotAdapterContent <- function(x, subset,
                               adapterType, plotType = "heatmap", nc,
                               pwfCols,
                               trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",  ...){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(trimNames))
  stopifnot(plotType %in% c("line", "heatmap"))

  # Get the AdapterContent
  x <- tryCatch(x[subset])
  df <- tryCatch(Adapter_Content(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Set the positions as a factor
  df$Position <- factor(df$Position, levels = unique(df$Position))

  # Change to long form and remove the _ symbols between words
  df <- reshape2::melt(df, id.vars = c("Filename", "Position"),
                       value.name = "Percent", variable.name = "Type")
  df <- dplyr::mutate(df, Type = gsub("_", " ", Type))

  # Restrict to a given adapterType if requested
  if (!missing(adapterType)) {
    keep <- grep(pattern = adapterType, x = df$Type, ...)
    df <- df[keep,]
    if(nrow(df) == 0) stop("No adapters matching the supplied type were found")
  }

  # Set the default number of columns for the faceting
  if (missing(nc)){
    # Choose 1 for the heatmap, 2 for the line
    nc <- grep(plotType, c("heatmap", "line"))
  }
  nc <- nc[1]
  stopifnot(is.numeric(nc))

  if (plotType == "heatmap"){

    acPlot <- ggplot(df, aes(x = Position, y = Filename, fill = Percent)) +
      geom_tile() +
      facet_wrap(~Type, ncol = nc)

    # Define the colour palette
    upr <- max(df$Percent)
    nCols <- findInterval(upr, c(0, 5, 10, Inf)) + 1
    gradCols <- getColours(pwfCols)[1:nCols]
    breaks <- c(0, 5, 10, upr)[1:nCols]

    # Add them to the initial plot
    acPlot <- acPlot +
      scale_fill_gradientn(colours = gradCols, values = breaks / max(breaks))
  }

  if (plotType == "line") {
    # Create the basic plot
    acPlot <- ggplot(df,
                     aes(x = as.integer(Position), y = Percent, colour = Filename)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 5,
               fill = getColours(pwfCols)["PASS"], alpha = 0.3) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5, ymax = 10,
               fill = getColours(pwfCols)["WARN"], alpha = 0.3) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = Inf,
               fill = getColours(pwfCols)["FAIL"], alpha = 0.3) +
      geom_line() +
      scale_y_continuous(limits = c(0, 100)) +
      scale_x_continuous(breaks = seq_along(levels(df$Position)),
                         labels = unique(df$Position),
                         expand = c(0, 0)) +
      ylab("Percent (%)")

    # Add the basic customisations
    acPlot <- acPlot + facet_wrap(~Type, ncol = nc[1])
  }

  # And draw the plot
  acPlot +
    labs(x = "Position in read (bp)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}
