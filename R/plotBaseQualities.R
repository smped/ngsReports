#' @title Plot the Base Qualities for each file
#'
#' @description Plot the Base Qualities for each file as separate plots
#'
#' @details This replicates the \code{Per base sequence quality} plots from FASTQC,
#' using facets to plce them all in a single ggplot2 object.
#'
#' For large datasets, subsetting by R1 or R2 reads may be helpful
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param nc \code{numeric}. The number of columns to create in the plot layout
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
#' # The default and subset plot
#' plotBaseQualities(fdl)
#' r1 <- grepl("R1", fileNames(fdl))
#' plotBaseQualities(fdl, subset = r1 )
#'
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate
#'
#' @export
plotBaseQualities <- function(x, subset, nc = 2,
                              trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  x <- x[subset]
  df <- tryCatch(Per_base_sequence_quality(x))
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start),
                      Start = as.factor(Start))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Set the y limit
  ylim <- c(0, max(df$`90th_Percentile`) + 1)

  # Get the Illumina encoding
  enc <- Basic_Statistics(x)$Encoding[1]
  enc <- gsub(".*(Illumina [0-9\\.]*)", "\\1", enc)

  qualPlot <- ggplot2::ggplot(df, ggplot2::aes(x = as.integer(Start), y = Median)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 30, ymax = Inf,
                      fill = rgb(0, 0.9, 0.6), alpha = 0.3) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30,
                      fill = rgb(0.9, 0.9, 0.7), alpha = 0.5) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 20,
                      fill = rgb(0.8, 0.4, 0.5), alpha = 0.5) +
    ggplot2::geom_crossbar(ggplot2::aes(ymin = Lower_Quartile, ymax = Upper_Quartile),
                           fill = "yellow", width = 0.8, size = 0.2) +
    ggplot2::geom_segment(ggplot2::aes(x = as.integer(Start)-0.4, xend = as.integer(Start) + 0.4,
                                       yend = Median), colour = "red") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = `10th_Percentile`, ymax = Lower_Quartile)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = Upper_Quartile, ymax = `90th_Percentile`)) +
    ggplot2::geom_line(ggplot2::aes(y = Mean), colour = "blue") +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(df$Start)),
                                labels = unique(df$Base),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0,0)) +
    ggplot2::xlab("Position in read (bp)") +
    ggplot2::ylab(paste0("Quality Scores (", enc, " encoding)")) +
    ggplot2::facet_wrap(~Filename, ncol = nc) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  # Draw the plot
  qualPlot

}
