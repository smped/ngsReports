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
#' @include validationFunctions.R
#'
#' @slot path Character vector of length 1 which contains a valid file path.
#' @export
setClass("FastqcFile", slots = c(path = "character"))
setValidity("FastqcFile", isValidFastqcFile)


#' @title Create a new FastqcFile Object
#' @description Create a new FastqcFile Object
#' @details Create a new FastqcFile Object from an external file
#' @return An object of class FastqcFile
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
