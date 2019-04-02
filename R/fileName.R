#' @title Return the Underlying Fastq File Names from FastqcData* Objects
#'
#' @description Return the Underlying Fastq File Names from FastqcData* Objects
#'
#' @param object An object of class FastqcData or FastqcDataList
#'
#' @return Returns the filenames for the Fastq file the FastQC report was
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
#' fileName(fdl)
#'

#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @name fileName
#' @aliases fileName,FastqcData-method
#' @rdname fileName-methods
setMethod(
    "fileName",
    "FastqcData",
    function(object){object@Summary$Filename[1]}
)

#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @name fileName
#' @aliases fileName,FastqcDataList-method
#' @rdname fileName-methods
setMethod(
    "fileName",
    "FastqcDataList",
    function(object){vapply(object@.Data, fileName, character(1))}
)

