#' @title Get the FASTQC version
#'
#' @description Get the FASTQC version used to generate the intial files
#'
#' @param object An object of class \code{FastqcData} or \code{FastqcDataList}
#'
#' @return A character vector (FastqcData), or tibble (FastqcDataList)
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Get the FASTQC version
#' fqcVersion(fdl)
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @docType methods

#' @export
#' @rdname fqcVersion
#' @aliases fqcVersion
setMethod("fqcVersion", "FastqcData", function(object){object@version})

#' @export
#' @rdname fqcVersion
#' @aliases fqcVersion
setMethod("fqcVersion", "FastqcDataList", function(object){
    tibble::tibble(
        Filename = fqName(object),
        version = vapply(object@.Data, fqcVersion, character(1))
    )
})

