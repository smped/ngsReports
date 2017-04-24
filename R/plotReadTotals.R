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
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param th A \code{ggplot2 theme} object. Defaults to \code{theme_bw()}
#'
#' @return Returns a ggplot object.
#'
#' @import ggplot2
#' @importFrom stringr str_detect
#'
#' @export
plotReadTotals <- function(x, subset, trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*", th){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  df <- readTotals(x, subset = subset, trimNames = trimNames, pattern = pattern)

  rtPlot <- ggplot2::ggplot(df, ggplot2::aes(x = Filename, y = Total_Sequences/1e06)) +
    ggplot2::geom_bar(stat = "identity")

  # Apply theme_bw() if missing, or an invalid theme is supplied
  if (missing(th)) {
    th <- ggplot2::theme_bw()
  }
  else{
    if (!ggplot2::is.theme(th)) {
      warning("Theme supplied as th is not a ggplot theme and will be ignored")
      th <- ggplot2::theme_bw()
    }
  }

  rtPlot +
    th +
    labs(y = "Total Reads (millions)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}
