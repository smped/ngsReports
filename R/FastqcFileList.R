#' @title The .FastqcFileList Object Class
#'
#' @description The .FastqcFileList Object Class
#'
#' @details The is an object which refers to a list of fastqc output files.
#' Only the paths are stored, however the files are checked for the correct
#' structure on formation of the object using \code{.FastqcFile}
#' Files can be any combination of zipped (*_fastqc.zip) or extracted
#' directories
#'
#' @param x character vector or generic list of .FastqcFiles
#' @return An object of class .FastqcFileList
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a .FastqcFileList
#' ffl <- .FastqcFileList(fileList)
#'
#' @include validationFunctions.R
#'
#' @slot ... this can either be a single character vector of paths to FASTQC
#' files, or several instances of .FastqcFile objects
#'
#' @importFrom methods new
#' @keywords internal
setClass(".FastqcFileList", contains = "list")
setValidity(".FastqcFileList", .isValidFastqcFileList)
.FastqcFileList <- function(x){

    if (is.character(x)) {
        fls <- lapply(x, .FastqcFile)
    }
    if (is.list(x)) {
        ## Check the class of every list element
        cls <- vapply(x, class, character(1))
        ## Exit if any elements are not .FastqcFiles
        if (any(!cls %in% ".FastqcFile")) {
            msg <-  paste0(
                "Method can only be applied to a generic list\n",
                "of .FastqcFile objects."
            )
            stop(msg)
        }
        fls <- x
    }

    new(".FastqcFileList", fls)
}


## The show method doesn't need exporting
setMethod("show", ".FastqcFileList",function(object){
    l <- length(object)
    cat(".FastqcFileList of", l, "file(s).\n")
    cat("Located in:\n",
        paste(unique(dirname(path(object))), collapse = "\n"))
})
