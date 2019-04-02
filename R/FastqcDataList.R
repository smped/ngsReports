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
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#' fdl
#'
#' @include validationFunctions.R
#'
#' @slot ... this can either be a single character vector of paths to FASTQC
#' files, or several instances of .FastqcFile objects
#' @rdname FastqcDataList
#' @aliases FastqcDataList-class
setClass("FastqcDataList", contains = "list")
setValidity("FastqcDataList", .isValidFastqcDataList)

#' @param x Character vector of file paths specifying paths to FastQC reports
#' @rdname FastqcDataList
#' @aliases FastqcDataList-class
#' @export
FastqcDataList <- function(x){

    fls <- lapply(x, .FastqcFile)
    fdl <- lapply(fls, as, "FastqcData")
    names(fdl) <- basename(x)
    new("FastqcDataList", fdl)
}

## The show method doesn't need exporting
setMethod(
    "show",
    "FastqcDataList",
    function(object){
        l <- length(object)
        cat("FastqcDataList for", l, "files.\n")
    }
)
