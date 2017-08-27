#' @title Draw a barplot of read totals
#'
#' @description Draw a barplot of read totals
#'
#' @details Draw a barplot of read totals using the standard ggplot2 syntax.
#' Read totals will be plotted in millions as this is the most common.
#' The raw data from \code{\link{readTotals}} can otherwise be used to manually create a plot.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param millions \code{logical}. Use Millions of reads as the scale for the y-axis.
#' Unless specified, will be set as TRUE automatically if the highest total is > 2e06.
#' @param ... Used to pass additional attributes to theme()
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
#' # Plot the Read Totals
#' plotReadTotals(fdl)
#'
#'
#' @return Returns a ggplot object.
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom plotly ggplotly
#'
#' @export
plotReadTotals <- function(x, usePlotly = FALSE, labels, millions, ...){

  df <- tryCatch(readTotals(x))

  # Automatically determine whether to convert to millions
  # Automatically determine whether to convert to millions
  if (missing(millions)) {
    millions <- ifelse(max(df$Total_Sequences) > 2e06, TRUE, FALSE)
  }
  millions <- millions[1]
  stopifnot(is.logical(millions))
  ylab <- c("Read Totals", "Read Totals (millions)")[millions + 1]

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)),
                        names = unique(df$Filename))
  }
  else{
    if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  # Setup the basic plot in millions or not
  if (millions){
    rtPlot <- ggplot(df, aes(x = Filename, y = Total_Sequences/1e06)) +
      labs(y = "Total Reads (millions)")
  }
  else{
    rtPlot <- ggplot(df, aes(x = Filename, y = Total_Sequences)) +
      labs(y = "Total Reads")
  }

  # Add the rest of the parameters
  rtPlot <- rtPlot +
    geom_bar(stat = "identity") +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(userTheme)) rtPlot <- rtPlot + userTheme

  if(usePlotly){
        rtPlot <- ggplotly(rtPlot)
  }
  # Draw the plot
  rtPlot

}
