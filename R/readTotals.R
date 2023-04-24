#' @title Get the read totals
#'
#' @description Get the read totals from one or more FASTQC reports
#'
#' @param x Can be a `FastqcData`, `FastqcDataList`, `FastpData`,
#' `FastpDataList` or file paths
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

  if (is(x, "FastqcDataList") | is(x, "fastqcData")) {
    df <-  tryCatch(getModule(x, "Basic_Statistics"))
    out <- df[c("Filename", "Total_Sequences")]
  }
  if (is(x, "FastpDataList")) {
    df <- tryCatch(getModule(x, "Summary"))$Before_filtering
    out <- df[c("Filename", "total_reads")]
    isPaired <- getModule(x, "paired")$paired
    out$total_reads[isPaired] <- out$total_reads[isPaired] / 2
  }
  if (is(x, "FastpData")) {
    df <- tryCatch(getModule(x, "Summary"))$Before_filtering
    out <- df[c("Filename", "total_reads")]
    isPaired <- getModule(x, "paired")
    if (isPaired) out$total_reads <- out$total_reads / 2
  }
  ## Tidy up any fastp column names
  names(out) <- gsub("total_reads", "Total_Sequences", names(out))
  out
}
