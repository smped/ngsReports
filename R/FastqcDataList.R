#' @title The FastqcDataList Object Class
#'
#' @description The FastqcDataList Object Class
#'
#' @return An object of class FastqcDataList
#'
#' @include validationFunctions.R
#'
#' @slot ... this can either be a single character vector of paths to FASTQC files, or several instances of FastqcFile objects
setClass("FastqcDataList", contains="list")
setValidity("FastqcDataList", isValidFastqcDataList) # Not written or defined yet

# The show method doesn't need exporting
setMethod("show", "FastqcDataList",
          function(object){
            l <- length(object)
            cat("FastqcDataList for", l, "files.\n")
          })
