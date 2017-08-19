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
#' @param warn,fail The default values for warn and fail are 5 and 10 respectively (i.e. precentages)
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
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
#' # Plot just the Universal Adapter
#' # and change the y-axis using ggplot2::scale_y_continuous
#' library(ggplot2)
#' plotAdapterContent(fdl, adapterType ="Universal", plotType = "line") +
#' scale_y_continuous()
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_colour_discrete
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#'
#' @export
plotAdapterContent <- function(x, subset,
                               adapterType, plotType = "heatmap",
                               warn = 5, fail = 10,
                               pwfCols,
                               labels){

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(plotType %in% c("line", "heatmap"))

  # Get the AdapterContent
  x <- tryCatch(x[subset])
  df <- tryCatch(Adapter_Content(x))

  # Sort out the colours & pass/warn/fail breaks
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))
  stopifnot(is.numeric(c(warn, fail)))
  stopifnot(all(fail < 100, warn < fail,  warn > 0))

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Set the positions as a factor
  df$Position <- factor(df$Position, levels = unique(df$Position))

  # Change to long form and remove the _ symbols between words
  df <- reshape2::melt(df, id.vars = c("Filename", "Position"),
                       value.name = "Percent", variable.name = "Type")
  df <- dplyr::mutate(df, Type = gsub("_", " ", Type))

  # Restrict to a given adapterType if requested
  if (!missing(adapterType)) {
    keep <- grep(pattern = adapterType, x = df$Type)
    df <- df[keep,]
    if(nrow(df) == 0) stop("No adapters matching the supplied type were found")
  }

  if (plotType == "heatmap"){

    # Define the colour palette
    upr <- max(df$Percent)
    nCols <- findInterval(upr, c(0, warn, fail, 100)) + 1
    gradCols <- getColours(pwfCols)[1:nCols]
    breaks <- c(0, warn, fail, 100)[1:nCols]

    # Reverse the factor levels for a better looking plot
    df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

    acPlot <- ggplot(df, aes_string(x = "Position", y = "Filename", fill = "Percent")) +
      scale_y_discrete(labels = labels) +
      geom_tile() +
      facet_wrap(~Type, ncol = 1) +
      scale_fill_gradientn(colours = gradCols, values = breaks / max(breaks))

  }

  if (plotType == "line") {
    df$x <- as.integer(df$Position)
    # Create the basic plot
    acPlot <- ggplot(df,
                     aes_string(x = "x", y = "Percent", colour = "Filename")) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = warn,
               fill = getColours(pwfCols)["PASS"], alpha = 0.3) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = warn, ymax = fail,
               fill = getColours(pwfCols)["WARN"], alpha = 0.3) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = fail, ymax = Inf,
               fill = getColours(pwfCols)["FAIL"], alpha = 0.3) +
      geom_line() +
      scale_y_continuous(limits = c(0, 100)) +
      scale_x_continuous(breaks = seq_along(levels(df$Position)),
                         labels = unique(df$Position),
                         expand = c(0, 0)) +
      scale_colour_discrete(labels = labels) +
      ylab("Percent (%)")

    # Add the basic customisations
    acPlot <- acPlot + facet_wrap(~Type, ncol = 1)
  }

  # And draw the plot
  acPlot +
    labs(x = "Position in read (bp)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}
