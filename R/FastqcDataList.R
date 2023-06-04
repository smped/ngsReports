#' @title The FastqcDataList Object Class
#'
#' @description The FastqcDataList Object Class
#'  `r lifecycle::badge("stable")`
#'
#' @slot ... this can either be a single character vector of paths to FASTQC
#' files, or several instances of .FastqcFile objects
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
#' @rdname FastqcDataList
#' @aliases FastqcDataList-class
setClass("FastqcDataList", contains = "list")
setValidity("FastqcDataList", .isValidFastqcDataList)

#' @param x Character vector of file paths specifying paths to FastQC reports
#' @rdname FastqcDataList
#' @aliases FastqcDataList-class
#' @export
FastqcDataList <- function(x){

    stopifnot(length(x) > 0)
    if (!any(is(x, "character"), is(x, "list"))) .errNotImp(x)

    if (is.character(x)) {
        fls <- lapply(x, .FastqcFile)
        fdl <- lapply(fls, as, "FastqcData")
        names(fdl) <- x
        fdl <- as(fdl, "FastqcDataList")
    }
    if (is.list(x)) fdl <- as(x, "FastqcDataList")
    fdl
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
