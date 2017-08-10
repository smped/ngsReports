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
setClass("FastqcFile", slots = c(path = "character"))
setValidity("FastqcFile", isValidFastqcFile)

#' Names cannot be set on this type of object. Not required for export
#' @param x An object of class FastqcFile
#' @param value Not used
setMethod("names<-", "FastqcFile", function(x, value){
  warning("The names attribute cannot be set on a FastqcFile object")
  x
})

# The show method doesn't need exporting
setMethod("show", "FastqcFile",
          function(object){
            cat(fileNames(object), "\n")
            cat("Located in", dirname(path(object)), "\n")
          })
