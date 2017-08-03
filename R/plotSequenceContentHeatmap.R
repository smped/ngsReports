#' @title Plot the per base quality as a heatmap
#'
#' @description Plot the Per Base Sequence Quality for a set of FASTQC files
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param type \code{character} Type of quality data to be presented "mean" or "median"
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param pwfcols Object of class \code{\link{Pwfcol}} to give colours for pass, warning, and fail
#' values in plot
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param
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
#' # The default plot
#' plotGcHeatmap(fdl)
#'
#' # Using counts
#' plotGcHeatmap(fdl, counts = TRUE)
#'
#' @import ggplot2
#' @import scales
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#'
#' @export

plotSequenceContent <- function(x, subset, pwfCols, clusterNames = TRUE){
  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  # Sort out the colours
  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  col <- ngsReports::getColours(pwfCols)

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }

  # Get the NContent
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_base_sequence_content(x))
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start))

  #get longest sequence
  basicStat <- Basic_Statistics(fdl) %>% dplyr::select(Filename, Longest_sequence)

  df <- df %>% dplyr::right_join(basicStat, by = "Filename") %>%
    dplyr::select(Filename, Start, , Longest_sequence) %>%
    tidyr::spread(Start, Mean)




}
















