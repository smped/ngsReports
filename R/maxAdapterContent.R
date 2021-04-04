#' @title Get the maximum Adapter Content
#'
#' @description Get the maximum Adapter Content across one or more FASTQC
#' reports
#'
#' @details This will extract the `Adapter_Content` module from the
#' supplied object, and provide a `tibble` with the final value for each
#' file.
#'
#' @param x Can be a `.FastqcFile`, `FastqcData`,
#' `FastqcDataList` or path
#' @param asPercent `logical`.
#' Format the values as percentages with the added `\%` symbol
#'
#' @return A `tibble` object containing the percent of reads with each
#' adapter type at the final position
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Get the maxAdapterContent
#' maxAdapterContent(fdl)
#'
#' @importFrom tidyselect one_of
#' @importFrom tidyr pivot_wider
#'
#' @export
maxAdapterContent <- function(x, asPercent = TRUE){

    stopifnot(is.logical(asPercent))

    ## Get the AdapterContent
    ac <- tryCatch(getModule(x, "Adapter_Content"))

    ## Perform the summary
    adapters <- setdiff(colnames(ac), c("Filename", "Position"))
    ac <- tidyr::gather(ac, "Type", "value", one_of(adapters))
    ac <- dplyr::group_by(ac, Filename, Type)
    ac <- dplyr::summarise_at(ac, dplyr::vars("value"), max)

    ## Format & return the output
    if (asPercent) ac$value <- scales::percent_format(0.01, 1)(ac$value)
    pivot_wider(ac, names_from = "Type", values_from = "value")

}
