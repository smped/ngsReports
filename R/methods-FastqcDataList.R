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
#' @rdname FastqcDataList-methods
#' @aliases path,FastqcDataList-methods
setMethod("path", "FastqcDataList", function(object){vapply(object@.Data, path, character(1))})

#' @rdname FastqcDataList-methods
#' @aliases names,FastqcDataList-methods
setMethod("names", "FastqcDataList", function(x){vapply(x@.Data, names, character(1))})

#' @export
#' @rdname FastqcDataList-methods
#' @aliases fileNames,FastqcDataList-method
setMethod("fileNames", "FastqcDataList", function(object){vapply(object@.Data, fileNames, character(1))})

#' @export
#' @rdname FastqcDataList-methods,
#' @aliases getSummary
setMethod("getSummary", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, getSummary)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname Version
setMethod("Version", "FastqcDataList",
          function(object){
            data_frame(Filename = fileNames(object),
                       Version = vapply(object@.Data, Version, character(1)))
          })

