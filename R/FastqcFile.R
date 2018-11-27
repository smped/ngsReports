#' @title The FastqcFile Object Class
#'
#' @description The FastqcFile Object Class
#'
#' @details The is an object which refers to a fastqc output file.
#' Only the path is stored, however the file is checked for the correct structure on formation of the object.
#' Files can be zipped (*_fastqc.zip) or extracted directories
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
#' # As this is the root structure, we can only call this
#' # function with an individual file
#' ff <- FastqcFile(fileList[[1]])
#'
#' @include validationFunctions.R
#'
#' @slot path Character vector of length 1 which contains a valid file path.
#' @export
#' @rdname FastqcFile
setClass("FastqcFile", slots = c(path = "character"))
setValidity("FastqcFile", isValidFastqcFile)

#' @param x Character vector (1) specifying a valid path to a file/directory as output by FastQC
#' @importFrom methods new
#' @export
#' @rdname FastqcFile
#' @aliases FastqcFile
setGeneric("FastqcFile",function(x){standardGeneric("FastqcFile")})

#' @export
#' @rdname FastqcFile
#' @aliases FastqcFile
setMethod("FastqcFile", "character", function(x){new("FastqcFile", path = x)})

# The show method doesn't need exporting
setMethod("show", "FastqcFile",
          function(object){
            cat(fileName(object), "\n")
            cat("Located in", dirname(path(object)), "\n")
          })
