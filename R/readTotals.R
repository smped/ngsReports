#' @title Get the read totals
#'
#' @description Get the read totals from one or more FASTQC reports
#'
#' @param x Can be a `FastqcData`, `FastqcDataList` or file paths
#'
#' @return A `tibble` with the columns `Filename` and
#' `Total_Sequences`
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Print the read totals
#' readTotals(fdl)
#'
#' @export
readTotals <- function(x){

    df <-  tryCatch(getModule(x, "Basic_Statistics"))
    df[c("Filename", "Total_Sequences")]

}
