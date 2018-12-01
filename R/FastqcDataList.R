#' @title The FastqcDataList Object Class
#'
#' @description The FastqcDataList Object Class
#'
#' @return An object of class FastqcDataList
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' fdl
#'
#' @include validationFunctions.R
#'
#' @slot ... this can either be a single character vector of paths to FASTQC
#' files, or several instances of FastqcFile objects
setClass("FastqcDataList", contains = "list")
setValidity("FastqcDataList", isValidFastqcDataList)

## The show method doesn't need exporting
setMethod(
    "show",
    "FastqcDataList",
    function(object){
        l <- length(object)
        cat("FastqcDataList for", l, "files.\n")
    }
)
