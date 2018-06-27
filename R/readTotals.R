#' @title Get the read totals
#'
#' @description Get the read totals from one or more FASTQC reports
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' 
#' @return A \code{data_frame} with the columns \code{Filename} and \code{Total_Sequences}
#' 
#' @examples 
#' 
#' # Get the files included with the package
#' fileDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(fileDir, pattern = "zip$", full.names = TRUE)
#' 
#' # Form a FastqcDataList
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the read totals
#' readTotals(fdl)    
#'
#' @export
readTotals <- function(x){

  df <-  tryCatch(Basic_Statistics(x))
  df[c("Filename", "Total_Sequences")]

}
