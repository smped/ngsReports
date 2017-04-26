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
#' @param adapterType A regular expression used to filter which adapter(s) are plotted
#' @param plotType \code{character}. Can only take the values \code{plotType = "heatmap"}
#' or \code{plotType = "line"}
#' @param nc Number of columns to use when faceting by Adapter Type.
#' @param pass The colour for percentages considered as a PASS
#' @param warn The colour for percentages consider as a WARN
#' @param fail The colour for percentages consider as a FAIL
#' @param trimNames \code{logical}. Remove the file suffix from the names displyed in the legend.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
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
#' fileList <- system.file("extdata", fileList, package = "fastqcReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot
#' plotAdapterContent(fdl)
#'
#' # Also subset the reads to just the R1 files
#' r1 <- grepl("R1", fileNames(fdl))
#' plotAdapterContent(fdl, subset = r1)
#'
#' # Plot just the Universal Adapter, and change the y-axis
#' plotAdapterContent(fdl, adapterType ="Universal", plotType = "line") +
#' scale_y_continuous()
#'
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom stringr str_detect
#'
#' @export
plotAdapterContent <- function(x, subset,
                               adapterType, plotType = "heatmap", nc,
                               pass = rgb(0, 0.8,0), warn = rgb(0.9, 0.9, 0.2),
                               fail = rgb(0.8, 0.2, 0.2),
                               trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*",  ...){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(plotType %in% c("line", "heatmap"))

  # Get the AdapterContent
  x <- x[subset]
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

    acPlot <- ggplot(df,
                     ggplot2::aes(x = Position, y = Filename, fill = Percent)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~Type, ncol = nc)

    # Define the colour palette
    upr <- max(df$Percent)
    nCols <- findInterval(upr, c(0, 5, 10, Inf)) + 1
    gradCols <- c(pass, warn, fail, rgb(1,1,1))[1:nCols]
    breaks <- c(0, 5, 10, upr)[1:nCols]

    # Add them to the initial plot
    acPlot <- acPlot +
      ggplot2::scale_fill_gradientn(colours = gradCols, values = breaks / max(breaks))
  }

  if (plotType == "line") {
    # Create the basic plot
    acPlot <- ggplot2::ggplot(df,
                              ggplot2::aes(x = as.integer(Position), y = Percent, colour = Filename)) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 5,
                        fill = pass, alpha = 0.3) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5, ymax = 10,
                        fill = warn, alpha = 0.3) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = Inf,
                        fill = fail, alpha = 0.3) +
      ggplot2::geom_line() +
      ggplot2::scale_y_continuous(limits = c(0, 100)) +
      ggplot2::scale_x_continuous(breaks = seq_along(levels(df$Position)),
                                  labels = unique(df$Position),
                                  expand = c(0, 0)) +
      ggplot2::ylab("Percent (%)")

    # Add the basic customisations
    acPlot <- acPlot + ggplot2::facet_wrap(~Type, ncol = nc[1])
  }

  # And draw the plot
  acPlot +
    ggplot2::labs(x = "Position in read (bp)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

}
