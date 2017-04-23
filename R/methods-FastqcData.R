#' Methods for an object of class FastqcData
#'
#' @description Methods and Accessors for an object of class FastqcData
#'
#' @details Describes all methods for an object of class \code{FastqcData}
#'
#' @param object An object of class \code{FastqcData}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @importFrom scales comma
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr everything
#'
#' @export
#' @rdname FastqcData-methods
#' @aliases path,FastqcData-method
setMethod("path", "FastqcData", function(object){object@path})

#' @export
#' @rdname FastqcData-methods
#' @aliases Version,FastqcData-method
setMethod("Version", "FastqcData", function(object){object@Version})

#' @export
#' @rdname FastqcData-methods
#' @aliases getSummary,FastqcData-method
setMethod("getSummary", "FastqcData", function(object){object@Summary})

#' @export
#' @rdname FastqcData-methods
#' @aliases names,FastqcData-method
setMethod("names", "FastqcData", function(x){basename(x@path)})

#' @export
#' @rdname FastqcData-methods
#' @aliases fileName,FastqcData-method
setMethod("fileNames", "FastqcData", function(object){object@Summary$FastqFile[1]})
