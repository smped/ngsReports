#' Get the FASTQC version
#'
#' @description Get the FASTQC version used to generate the intial files
#'
#'
#' @param object An object of class \code{FastqcData} or \code{FastqcDataList}
#'
#' @return A character vector
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods

#' @export
#' @rdname Version
#' @aliases Version
setMethod("Version", "FastqcData", function(object){object@Version})

#' @export
#' @rdname Version
#' @aliases Version
setMethod("Version", "FastqcDataList",
          function(object){
            data_frame(Filename = fileName(object),
                       Version = vapply(object@.Data, Version, character(1)))
          })

