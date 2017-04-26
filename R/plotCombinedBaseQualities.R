#' @title Plot the combined Per_base_sequence_quality information
#'
#' @description Plot the Per_base_sequence_quality information for a set of FASTQC reports
#'
#' @details This enables plotting of any of the supplied values (i.e. Mean, Median etc) across
#' a set of FASTQC reports.
#' By default only the Mean will be plotted,
#' however any number of the supplied values can be added to the plot,
#' and these will be separated by linetype.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param value \code{character}. Specify which value whould be plotted.
#' Can be any of the columns returned by \code{\link{Per_base_sequence_quality}}.
#' Defaults to \code{value = "Mean"}.
#' Can additionally set to "all" to plot all available quantities
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
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
#' # Find the R1 files
#' r1 <- grepl("R1", fileNames(fdl))
#'
#' # The default plot using the Mean only
#' plotCombinedBaseQualities(fdl)
#'
#' # Plot the R1 files showing the Mean and Lower_Quartile
#' plotCombinedBaseQualities(fdl, subset = r1, value = c("Mean", "Lower_Quartile"))
#'
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom reshape2 melt
#'
#' @export
plotCombinedBaseQualities <- function(x, subset, value = "Mean",
                                      trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(is.character(value))

  x <- x[subset]
  df <- tryCatch(Per_base_sequence_quality(x))

  # Check for valid columns
  value <- setdiff(value, c("Filename", "Base")) # Not relevant columns
  if (length(value) == 0) {
    message("Invalid column specified, setting as 'Mean' (default)")
    value <- "Mean"
  }
  # Any vector containing 'all' will return all values...
  if (!any(grepl("all", value))){
    value <- intersect(value, colnames(df)[-c(1:2)])
    if (length(value) == 0) stop("The specified value could not be found in the output from Per_base_sequence_quality")
  }

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the Illumina encoding
  enc <-  Basic_Statistics(x)$Encoding[1]
  enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)

  # Find the central position for each base as some may be grouped
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start),
                      End = gsub("[0-9]*-([0-9]*)", "\\1", Base),
                      End = as.integer(End),
                      Base = 0.5*(Start + End))
  df <- dplyr::select(df, -Start, -End)
  if (!any(grepl("all", value))) df <- dplyr::select(df, Filename, Base, dplyr::one_of(value))
  df <- reshape2::melt(df, id.vars = c("Filename", "Base"),
                       variable.name = "Value", value.name = "Score")

  # Make basic plot, adding the shaded background colours
  qualPlot <- ggplot2::ggplot(df, ggplot2::aes(x = Base, y = Score, colour = Filename)) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = 30, ymax = Inf,
                      fill = rgb(0, 0.9, 0.6), alpha = 0.3) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = 20, ymax = 30,
                      fill = rgb(0.9, 0.9, 0.7), alpha = 0.5) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = 20,
                      fill = rgb(0.8, 0.4, 0.5), alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(linetype = Value)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::ylab(paste0("Quality Scores (", enc, " encoding)")) +
    ggplot2::theme_bw()

  # Draw the plot
  qualPlot

}
