#' Methods for an object of class StarLogData
#'
#' @description Methods and Accessors for an object of class StarLogData
#'
#' @details Describes all methods for an object of class \code{StarLogData}
#'
#' @param object An object of class \code{StarLogData}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @export
setMethod("fileNames", "StarLogData", function(object){object@fileName})

setMethod("show", "StarLogData",
          function(object){
            n <- length(fileNames(object))
            cat("Log data for", n, "files after STAR alignment")
          })
