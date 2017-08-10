#' @title Return the File Paths from an object
#'
#' @description Return the File Paths from an object
#'
#' @details Obtains the file.path for objects of multiple classes
#'
#' @param object An object of class FastqcFile or FastqcFileList
#'
#' @return A character(1) vector of the file name.
#'
#' @importMethodsFrom Rsamtools path
#' @export
setMethod("path", "FastqcFile", function(object){object@path})

#' @export
setMethod("path", "FastqcFileList", function(object){vapply(object, path, character(1))})

#' @export
setMethod("path", "FastqcData", function(object){object@path})

#' @export
setMethod("path", "FastqcDataList", function(object){vapply(object@.Data, path, character(1))})
