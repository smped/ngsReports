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
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileName
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName when \code{trimNames = TRUE}.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @importFrom stringr str_detect
#'
#' @export
readTotals <- function(x, subset, trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

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

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Return the data
  df

}
