#' @title The FastqcFile Object Class
#'
#' @description The FastqcFile Object Class defines a path to the output from
#' the standalone tool FastQC
#'
#' @details The is an object which refers to a fastqc output file.
#' Only the path is stored, however the file is checked for the correct
#' structure on formation of the object.
#' Files can be zipped (*_fastqc.zip) or extracted directories
#'
#'
#' @return An object of class FastqcFile
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # As this is the root structure, we can only call this
#' # function with an individual file
#' ff <- FastqcFile(fileList[[1]])
#'
#' @include validationFunctions.R
#'
#' @seealso \code{\link{FastqcFileList}}
#'
#' @slot path Character vector of length 1 which contains a valid file path.
#' @export
#' @rdname FastqcFile
setClass("FastqcFile", slots = c(path = "character"))
setValidity("FastqcFile", .isValidFastqcFile)

#' @param x character(1) denoting a file.path
#' @export
#' @rdname FastqcFile
#' @aliases FastqcFile
#' @importFrom methods new
FastqcFile <- function(x){
    ## Ensure only a single file/directory that exists is parsed
    stopifnot(is.character(x), length(x) == 1)
    stopifnot(file.exists(x))
    new("FastqcFile", path = x)
}

## The show method doesn't need exporting
setMethod("show", "FastqcFile", function(object){
    p <- path(object)
    cat(basename(p), "\n")
    cat("Located in", dirname(p), "\n")
})

