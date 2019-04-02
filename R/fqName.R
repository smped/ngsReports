#' @title Return the Underlying Fastq File Names from FastqcData* Objects
#'
#' @description Return the Underlying Fastq File Names from FastqcData* Objects
#'
#' @param object An object of class FastqcData or FastqcDataList
#'
#' @return Returns the names of the Fastq files the FastQC report was
#' generated from, without any preceding directories.
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#' fqName(fdl)
#'

#' @export
#' @name fqName
#' @aliases fqName,FastqcData-method
#' @rdname fqName-methods
setMethod(
    "fqName",
    "FastqcData",
    function(object){object@Summary$Filename[1]}
)

#' @export
#' @name fqName
#' @aliases fqName,FastqcDataList-method
#' @rdname fqName-methods
setMethod(
    "fqName",
    "FastqcDataList",
    function(object){vapply(object@.Data, fqName, character(1))}
)

