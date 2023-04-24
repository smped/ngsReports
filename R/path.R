#' @title Return the File Paths from an object
#'
#' @description Return the File Paths from an object
#'
#' @details Obtains the file.path for objects of multiple classes
#'
#' @param object An object of class .FastqcFile
#'
#' @return A character vector of the file paths to the underlying FastQC reports
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#' path(fdl)
#'
#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,.FastqcFile-method
#' @export
setMethod("path", ".FastqcFile", function(object){object@path})

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

#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,.FastpFile-method
#' @export
setMethod("path", ".FastpFile", function(object){object@path})

#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,FastpData-method
#' @export
setMethod("path", "FastpData", function(object){object@path})

#' @importMethodsFrom BiocGenerics path
#' @name path
#' @aliases path,FastpDataList-method
#' @export
setMethod("path", "FastpDataList", function(object){
  vapply(object@.Data, path, character(1))
})
