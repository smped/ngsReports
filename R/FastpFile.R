#' @title The .FastpFile Object Class
#'
#' @description The .FastpFile Object Class defines a path to the output from
#' the standalone tool fastp.
#'  `r lifecycle::badge("experimental")`
#'
#' @details This class simply refers to a fastp output file after checking for
#' existence and validity (i.e. the correct internal structure).
#' Underlying files are expected to be in json format
#'
#' The helper function `.FastpFile()` is a simple constructor which
#' checks validity and enables construction of other dependent classes.
#'
#' @param x character(1) denoting a file.path
#'
#' @slot path Character vector of length 1 which contains a valid file path.
#'
#' @return An object of class .FastqcFile
#'
#' @include validationFunctions.R
#'
#' @importFrom methods new
#' @keywords internal
setClass(".FastpFile", slots = c(path = "character"))
setValidity(".FastpFile", .isValidFastpFile)
.FastpFile <- function(x){

    ## Ensure only a single file/directory that exists is parsed
    stopifnot(!is.null(x))
    stopifnot(is.character(x), length(x) == 1)
    stopifnot(file.exists(x))
    new(".FastpFile", path = x)

}

## The show method doesn't need exporting
setMethod("show", ".FastpFile", function(object){
    p <- path(object)
    cat(basename(p), "\n")
    cat("Located in", dirname(p), "\n")
})

