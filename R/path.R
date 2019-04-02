#' @title Return the File Paths from an object
#'
#' @description Return the File Paths from an object
#'
#' @details Obtains the file.path for objects of multiple classes
#'
#' @param object An object of class .FastqcFile or .FastqcFileList
#'
#' @return A character vector of the file paths to the underlying FastQC reports
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' path(fdl)
#'
#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,.FastqcFile-method
#' @export
setMethod("path", ".FastqcFile", function(object){object@path})

#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,.FastqcFileList-method
#' @export
setMethod("path", ".FastqcFileList", function(object){
    vapply(object, path, character(1))
})

#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,FastqcData-method
#' @export
setMethod("path", "FastqcData", function(object){object@path})

#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,FastqcDataList-method
#' @export
setMethod("path", "FastqcDataList", function(object){
    vapply(object@.Data, path, character(1))
})
