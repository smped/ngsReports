#' @title Extract Elements
#'
#' @description Extract elements from FastqcDataList Object
#'
#' @details Extract elements in a consistent manner with R conventions
#'
#' @param x A FastqcDataList or FastpDataList
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
#' @docType methods
#' @aliases [,FastqcDataList,numeric,missing-method
#' @rdname extract-methods
#' @export
setMethod(
    f = "[",
    signature = c(x = "FastqcDataList", i = "numeric", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        x <- x@.Data[as.integer(i)]
        if (length(x) == 0) stop("Object cannot be of length 0")
        new("FastqcDataList", x)
    }
)
#' @aliases [,FastqcDataList,character,missing-method
#' @rdname extract-methods
setMethod(
    f = "[",
    signature = c(x = "FastqcDataList", i = "character", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        i <- match(i, names(x))
        if (anyNA(i)) stop("One or more supplied names not found in object")
        x[i]
    }
)
#' @aliases [,FastqcDataList,logical,missing-method
#' @rdname extract-methods
setMethod(
    f = "[",
    signature = c(x = "FastqcDataList", i = "logical", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        if (length(i) != length(x)) stop("x & i are of different lengths")
        x <- x@.Data[i]
        if (length(x) == 0) stop("Object cannot be of length 0")
        new("FastqcDataList", x)
    }
)
#' @aliases [,FastqcDataList,ANY,missing-method
#' @rdname extract-methods
setMethod(
    f = "[",
    signature = c(x = "FastqcDataList", i = "ANY", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        cl <- class(i)
        message("Subsetting not implemented for objects of class ", cl)
    }
)
#' @aliases [,FastpDataList,numeric,missing-method
#' @rdname extract-methods
#' @export
setMethod(
    f = "[",
    signature = c(x = "FastpDataList", i = "numeric", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        x <- x@.Data[as.integer(i)]
        if (length(x) == 0) stop("Object cannot be of length 0")
        new("FastpDataList", x)
    }
)
#' @aliases [,FastpDataList,character,missing-method
#' @rdname extract-methods
setMethod(
    f = "[",
    signature = c(x = "FastpDataList", i = "character", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        i <- match(i, names(x))
        if (anyNA(i)) stop("One or more supplied names not found in object")
        x[i]
    }
)
#' @aliases [,FastpDataList,logical,missing-method
#' @rdname extract-methods
setMethod(
    f = "[",
    signature = c(x = "FastpDataList", i = "logical", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        if (length(i) != length(x)) stop("x & i are of different lengths")
        x <- x@.Data[i]
        if (length(x) == 0) stop("Object cannot be of length 0")
        new("FastpDataList", x)
    }
)
#' @aliases [,FastpDataList,ANY,missing-method
#' @rdname extract-methods
setMethod(
    f = "[",
    signature = c(x = "FastpDataList", i = "ANY", j = "missing"),
    definition = function(x, i, j, ..., drop = TRUE){
        cl <- class(i)
        message("Subsetting not implemented for objects of class ", cl)
    }
)

