#' @title Get the read totals
#'
#' @description Get the read totals from one or more FASTQC reports
#'
#' @return A \code{data_frame} with the columns \code{Filename} and \code{Total_Sequences}
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful for only return totals from R1 files, or any other subset
#' @param trimNames \code{logical}. Remove the suffix from fileNames
#'
#' @export
readTotals <- function(x, subset, trimNames = TRUE){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))

  # Get the data
  x <- x[subset]
  df <-  tryCatch(Basic_Statistics(x))
  df <- df[c("Filename", "Total_Sequences")]
  if (trimNames)  df$Filename <- gsub("\\.(fastq|fq).*", "", df$Filename)

  # Return the data
  df

}
