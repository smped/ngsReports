#' @title Methods for an object of class FastqcDataList
#'
#' @description Methods and Accessors for an object of class FastqcDataList
#'
#' @details Describes all methods for an object of class \code{FastqcDataList}
#'
#' @param object An object of class \code{FastqcDataList}
#'
#' @return Returns a \code{data_frame} for all the methods except \code{path}, \code{names},
#' \code{fileNames}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @importFrom dplyr bind_rows
#'
#' @export
setMethod("path", "FastqcDataList", function(object){vapply(object@.Data, path, character(1))})

#' @export
setMethod("fileNames", "FastqcDataList", function(object){vapply(object@.Data, fileNames, character(1))})

#' @export
setMethod("getSummary", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, getSummary)
            dplyr::bind_rows(df)
          })

#' @export
setMethod("Version", "FastqcDataList",
          function(object){
            data_frame(Filename = fileNames(object),
                       Version = vapply(object@.Data, Version, character(1)))
          })

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
