#' Get all data from FastQC files
#'
#' @description Read the information from the \code{fastqc_data.txt} files in
#' each .FastqcFile
#'
#' @param object Can be a .FastqcFile or .FastqcFileList, or paths to files
#'
#' @return An object of \code{FastqcData} or a \code{FastqcDataList}
#'
#' @include AllGenerics.R
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' @export
#' @rdname getFastqcData-methods
setGeneric("getFastqcData", function(object){standardGeneric("getFastqcData")})

#' @importFrom utils unzip
#' @name getFastqcData
#' @aliases getFastqcData,.FastqcFile-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", ".FastqcFile", function(object){

    as(object, "FastqcData")

})

#' @name getFastqcData
#' @aliases getFastqcData,NULL-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "NULL", function(object){
    if (is.null(object)) stop("No files have been provided.")
}
)
#' @name getFastqcData
#' @aliases getFastqcData,character-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "character", function(object){
    x <- .FastqcFileList(object)
    n <- length(x)
    if (n == 1) x[[1]]
    getFastqcData(x)
}
)
#' @name getFastqcData
#' @aliases getFastqcData,.FastqcFileList-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", ".FastqcFileList", function(object){
    as(object, "FastqcDataList")
}
)

