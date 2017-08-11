#' @title The FastqcFileList Object Class
#'
#' @description The FastqcFileList Object Class
#'
#' @details The is an object which refers to a list of fastqc output file.
#' Only the paths are stored, however the files are checked for the correct structure on formation of the object.
#' Files can be any combination of zipped (*_fastqc.zip) or extracted directories
#'
#'
#' @return An object of class FastqcFile
#'
#' @include validationFunctions.R
#'
#' @slot ... this can either be a single character vector of paths to FASTQC files, or several instances of FastqcFile objects
setClass("FastqcFileList", contains="list")
setValidity("FastqcFileList", isValidFastqcFileList)

#' @title Create a new FastqcFileList Object
#' @description Create a new FastqcFileList Object
#' @details Create a new FastqcFileList Object from a vector of file paths
#' @return An object of class FastqcFileList
#' @param x Character vector specifying a valid paths to files/directories as output by FastQC
#' @importFrom methods new
#' @export
#' @rdname FastqcFileList
#' @aliases FastqcFileList
setGeneric("FastqcFileList", function(path){standardGeneric("FastqcFileList")})
setMethod("FastqcFileList", "character",
          function(path)
          {
            fls <- lapply(path, FastqcFile)
            names(fls) <- basename(path)
            new("FastqcFileList", fls)
          })

#' @export
#' @rdname FastqcFileList
#' @aliases FastqcFileList
setMethod("FastqcFileList", "list",
          function(path)
          {
            cls <- vapply(path, class, character(1))
            if (any(!cls %in% "FastqcFile")) stop("Method can only be applied to\nFastqcFile objects as a generic list.")
            new("FastqcFileList", path)
          })

#' @export
setMethod("[", "FastqcFileList", function(x, i, j, ..., drop = TRUE){FastqcFileList(x@.Data[i])})

#' @export
setMethod("[[", "FastqcFileList", function(x, i, j){x@.Data[[i]]})

# The show method doesn't need exporting
setMethod("show", "FastqcFileList",
          function(object){
            l <- length(object)
            cat("FastqcFileList of", l, "file(s).\n")
            cat("Located in:\n", paste(unique(dirname(path(object))), collapse = "\n"))
          })
