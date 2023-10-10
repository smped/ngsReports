#' @title The FastpDataList Object Class
#'
#' @description The FastpDataList Object Class
#'  `r lifecycle::badge("stable")`
#'
#' @slot ... this can either be a single character vector of paths to fastp
#' files, or several instances of .FastpFile objects
#'
#' @return An object of class FastpDataList
#'
#' @include validationFunctions.R
#'
#' @rdname FastpDataList
#' @aliases FastpDataList-class
setClass("FastpDataList", contains = "list")
setValidity("FastpDataList", .isValidFastpDataList)

#' @param x Character vector of file paths specifying paths to fastp.json.gz output
#' @rdname FastpDataList
#' @aliases FastpDataList-class
#' @export
FastpDataList <- function(x){

    stopifnot(length(x) > 0)
    if (!any(is(x, "character"), is(x, "list"))) .errNotImp(x)

    if (is.character(x)) {
        fls <- lapply(x, .FastpFile)
        fdl <- lapply(fls, as, "FastpData")
        names(fdl) <- x
        fdl <- as(fdl, "FastpDataList")
    }
    if (is.list(x)) fdl <- as(x, "FastpDataList")
    fdl
}

## The show method doesn't need exporting
setMethod(
    "show",
    "FastpDataList",
    function(object){
        l <- length(object)
        cat("FastpDataList for", l, "files.\n")
    }
)
