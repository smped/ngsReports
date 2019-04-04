#' @title The .FastqcFile Object Class
#'
#' @description The .FastqcFile Object Class defines a path to the output from
#' the standalone tool FastQC.
#'
#' The helper function \code{.FastqcFile()} is a simple constructor which
#' checks validity and enables construction of other dependent classes.
#'
#' @details This class simply refers to a fastqc output file after checking for
#' existence and validity (i.e. the correct internal structure).
#' Underlying files can be zipped (*_fastqc.zip) or extracted directories
#'
#' @param x character(1) denoting a file.path
#'
#' @slot path Character vector of length 1 which contains a valid file path.
#'
#' @return An object of class .FastqcFile
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1]
#'
#' # As this is the root structure, we can only call this
#' # function with an individual file
#' ff <- ngsReports:::.FastqcFile(fl)
#'
#' @include validationFunctions.R
#'
#' @importFrom methods new
#' @keywords internal
setClass(".FastqcFile", slots = c(path = "character"))
setValidity(".FastqcFile", .isValidFastqcFile)
.FastqcFile <- function(x){

    ## Ensure only a single file/directory that exists is parsed
    stopifnot(!is.null(x))
    stopifnot(is.character(x), length(x) == 1)
    stopifnot(file.exists(x))
    new(".FastqcFile", path = x)
}

## The show method doesn't need exporting
setMethod("show", ".FastqcFile", function(object){
    p <- path(object)
    cat(basename(p), "\n")
    cat("Located in", dirname(p), "\n")
})

