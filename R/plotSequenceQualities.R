#' @title Plot the Per Sequence Quality Scores
#'
#' @description Plot the Per Sequence Quality Scores for a set of FASTQC reports
#'
#' @details Plots the distribution of average sequence quality scores across the set of files.
#' Values can be plotted either as counts (\code{counts = TRUE}) or as frequencies (\code{counts = FALSE}).
#'
#' Any faceting or scale adjustment can be performed after generation of the initial plot,
#' using the standard methods of ggplot2 as desired.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours for PASS/WARN/FAIL
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @return A ggplot2 object
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
#' # Draw the defualt plot
#' plotSequenceQualities(fdl)
#'
#' # Get the R1 files
#' r1 <- grepl("R1", fileNames(fdl))
#'
#' # Customise using the R1 subset, plotting counts,
#' # and faceting after the initial function call
#' plotSequenceQualities(fdl, subset = r1, counts = TRUE) +
#'   facet_wrap(~Filename, ncol = 2)
#'
#' 
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#'
#' @export
plotSequenceQualities <- function(x, subset, counts = FALSE, pwfCols,
                                  trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  stopifnot(grepl("(Fastqc|character)", class(x)))
  stopifnot(is.logical(trimNames))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  stopifnot(isValidPwf(pwfCols))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_sequence_quality_scores(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Get the correct y-axis label
  ylab <- c("Frequency", "Count")[counts + 1]

  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>%
      dplyr::mutate(Freq = Count / sum(Count)) %>%
      dplyr::ungroup()
    qualPlot <- ggplot(df, aes(x = Quality, y = Freq, colour = Filename))

  }
  else{

    # Initialise the plot using counts
    qualPlot <- ggplot(df, aes(x = Quality, y = Count, colour = Filename))

  }

  qualPlot <- qualPlot +
    annotate("rect", xmin = 30, xmax = Inf, ymin = -Inf, ymax = Inf,
                      fill = getColours(pwfCols)["PASS"], alpha = 0.3) +
    annotate("rect", xmin = 20, xmax = 30, ymin = -Inf, ymax = Inf,
                      fill = getColours(pwfCols)["WARN"], alpha = 0.3) +
    annotate("rect", xmin = -Inf, xmax = 20,, ymin = -Inf, ymax = Inf,
                      fill = getColours(pwfCols)["FAIL"], alpha = 0.3) +
    scale_x_continuous(limits = c(0, 41), expand = c(0, 0), breaks = seq(0, 40, by = 10)) +
    geom_line() +
    xlab("Mean Sequence Quality Per Read (Phred Score)") +
    ylab(ylab) +
    theme_bw()

  # Draw the plot
  qualPlot

}
