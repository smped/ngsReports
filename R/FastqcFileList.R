#' @title The FastqcFileList Object Class
#'
#' @description The FastqcFileList Object Class
#'
#' @details The is an object which refers to a list of fastqc output files.
#' Only the paths are stored, however the files are checked for the correct
#' structure on formation of the object.
#' Files can be any combination of zipped (*_fastqc.zip) or extracted
#' directories
#'
#'
#' @return An object of class FastqcFile
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcFileList
#' ffl <- FastqcFileList(fileList)
#'
#' @include validationFunctions.R
#'
#' @rdname FastqcFileList
#' @slot ... this can either be a single character vector of paths to FASTQC
#' files, or several instances of FastqcFile objects
setClass("FastqcFileList", contains = "list")
setValidity("FastqcFileList", isValidFastqcFileList)

#' @param x Character vector specifying a valid paths to files/directories as
#' output by FastQC
#' @importFrom methods new
#' @export
#' @rdname FastqcFileList
#' @aliases FastqcFileList
setGeneric("FastqcFileList", function(x){standardGeneric("FastqcFileList")})

#' @export
#' @rdname FastqcFileList
#' @aliases FastqcFileList
setMethod(
    "FastqcFileList",
    "character",
    function(x) {
        fls <- lapply(x, FastqcFile)
        new("FastqcFileList", fls)
    }
)

#' @export
#' @rdname FastqcFileList
#' @aliases FastqcFileList
setMethod(
    "FastqcFileList",
    "list",
    function(x) {
        cls <- vapply(x, class, character(1))
        if (any(!cls %in% "FastqcFile")) {
            msg <- paste0("Method can only be applied to\n",
                          "FastqcFile objects as a generic list.")
            stop(msg)
        }
        new("FastqcFileList", x)
    }
)

## The show method doesn't need exporting
setMethod(
    "show",
    "FastqcFileList",
    function(object){
        l <- length(object)
        cat("FastqcFileList of", l, "file(s).\n")
        cat("Located in:\n",
            paste(unique(dirname(path(object))), collapse = "\n"))
    }
)
