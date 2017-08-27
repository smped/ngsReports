#' @title Get the read totals
#'
#' @description Get the read totals from one or more FASTQC reports
#'
#' @return A \code{data_frame} with the columns \code{Filename} and \code{Total_Sequences}
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#'
#' @export
readTotals <- function(x){

  df <-  tryCatch(Basic_Statistics(x))
  df[c("Filename", "Total_Sequences")]

}
