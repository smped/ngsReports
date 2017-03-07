#' Define object classes
#'
FastqcFile <- setClass("FastqcFile", slots = c(path = "character", compressed = "logical"))
FastqcFileList <- setClass("FastqcFileList", contains="list")

#' @include validationFunctions.R
setValidity("FastqcFile", validFastqcFile)

# Set the Generics
setGeneric("FastqcFile",function(filePath){standardGeneric("FastqcFile")})
setGeneric("FastqcFileList", function(path){standardGeneric("FastqcFileList")})
setGeneric("isCompressed", function(object){standardGeneric("isCompressed")})
setGeneric("path", function(object){standardGeneric("path")})
