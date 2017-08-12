#' @title Methods for an object of class FastqcDataList
#'
#' @description Methods and Accessors for an object of class FastqcDataList
#'
#' @details Describes all methods for an object of class \code{FastqcDataList}
#'
#' @param object An object of class \code{FastqcDataList}
#'
#' @return Returns a \code{data_frame} for all the methods except \code{path}, \code{names},
#' \code{fileName}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @importFrom dplyr data_frame
#'

# Define the subsetting method, to remain consistent with list behaviour
setMethod("[", "FastqcDataList",
          function(x, i, j, ..., drop = TRUE){
            new("FastqcDataList", x@.Data[i])
          })

# The show method doesn't need exporting
setMethod("show", "FastqcDataList",
          function(object){
            l <- length(object)
            cat("FastqcDataList for", l, "file(s).\n")
          })
