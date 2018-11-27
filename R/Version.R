#' Get the FASTQC version
#'
#' @description Get the FASTQC version used to generate the intial files
#'
#'
#' @param object An object of class \code{FastqcData} or \code{FastqcDataList}
#'
#' @return A character vector (FastqcData), or data_frame (FastqcDataList)
#' 
#' @examples 
#' 
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' 
#' # Get the FASTQC version
#' Version(fdl)
#'
#' @include FastqcData.R
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
            tibble::tibble(Filename = fileName(object),
                           Version = vapply(object@.Data, Version, character(1)))
          })

