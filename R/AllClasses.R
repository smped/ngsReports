# Define object classes
FastqcFile <- setClass("FastqcFile", slots = c(path = "character", compressed = "logical"))
FastqcFileList <- setClass("FastqcFileList", contains="list")

# Set the vailation functions for any object classes
#' @include validationFunctions.R
setValidity("FastqcFile", validFastqcFile)

# Set the Generics
setGeneric("FastqcFile",function(filePath){standardGeneric("FastqcFile")})
setGeneric("FastqcFileList", function(path){standardGeneric("FastqcFileList")})
setGeneric("isCompressed", function(object){standardGeneric("isCompressed")})
setGeneric("path", function(object){standardGeneric("path")})
setGeneric("getSummary", function(object){standardGeneric("getSummary")})
setGeneric("getFastqcData", function(object){standardGeneric("getFastqcData")})
