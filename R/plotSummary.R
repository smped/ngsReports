#' @title Plot the PASS/WARN/FAIL information
#'
#' @description Extract the PASS/WARN/FAIL summaries and plot them
#'
#' @details This uses the standard ggplot2 syntax to create a three colour plot.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful for only return totals from R1 files, or any other subset
#' @param col \code{character vector} of colours
#' @param trimNames \code{logical}. Remove the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @return A standard ggplot2 object
#'
#' @import ggplot2
#'
#' @export
plotSummary <- function(x, subset, col = c(FAIL="red", WARN = "yellow", PASS="green"),
                        trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  x <- x[subset]
  df <- tryCatch(getSummary(x))
  if (trimNames) {
    stopifnot(is.character(pattern))
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
  }
  df$Category <- factor(df$Category, levels = rev(unique(df$Category)))
  df$Status <- factor(df$Status, levels = c("PASS", "WARN", "FAIL"))

  ggplot2::ggplot(df, aes(x = Filename, y = Category, fill = Status)) +
    geom_tile(colour = "black") +
    scale_fill_manual(values = col) +
    labs(x="Filename", y="QC Category") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}
