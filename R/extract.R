#' @title Extract Elements
#'
#' @description Extract elements from FastqcDataList Object
#'
#' @details Extract elements in a consistent manner with R conventions
#'
#' @param x A FastqcDataList
#' @param i character, logical or integer vector
#' @param j not used
#' @param ... not used
#' @param drop not used
#'
#' @include AllGenerics.R
#'
#' @return
#' Will return a subset of the original object following the standard
#' rules for subsetting objects
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
#' # Subsetting using the standard methods
#' fdl[1]
#' fdl[[1]]
#'
#' @name [
#' @aliases [,FastqcDataList,ANY,missing,ANY-method
#' @rdname extract-methods
#' @export
setMethod(
    f = "[",
    signature = c(x = "FastqcDataList", i = "ANY", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        new("FastqcDataList", x@.Data[i])
    }
)

