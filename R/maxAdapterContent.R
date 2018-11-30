#' @title Get the maximum Adapter Content
#'
#' @description Get the maximum Adapter Content across one or more FASTQC 
#' reports
#'
#' @details This will extract the \code{Adapter_Content} from the supplied 
#' object, and provide a \code{tibble} with the final value for each file
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, 
#' \code{FastqcData}, \code{FastqcDataList} or path
#' @param asPercent \code{logical}.
#' Format the values as percentages with the added \code{\%} symbol
#'
#' @return A \code{tibble} object containing the percent of reads with each 
#' adapter type at the final position
#' 
#' @examples 
#' # Get the files included with the package 
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#' 
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' 
#' # Get the maxAdapterContent
#' maxAdapterContent(fdl)
#'
#'
#' @export
maxAdapterContent <- function(x, asPercent = TRUE){
    
    stopifnot(is.logical(asPercent))
    
    # Get the AdapterContent
    ac <- tryCatch(Adapter_Content(x))
    
    # Perform the summary
    adapters <- setdiff(colnames(ac), c("Filename", "Position"))
    ac <- tidyr::gather(ac, "Type", "value", tidyselect::one_of(adapters))
    ac <- dplyr::group_by(ac, Filename, Type)
    ac <- dplyr::summarise_at(ac, dplyr::vars("value"), dplyr::funs("max"))
    
    # Format & return the output
    if (asPercent) ac$value <- scales::percent_format(0.01, 1)(ac$value)
    tidyr::spread(ac, "Type", "value")
    
}
