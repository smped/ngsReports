#' Get the set of Basic Statistics
#'
#' @description Retrieve the Basic Statistics module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#' @include FastqcFile.R
#' @include FastqcFileList.R
#' @include FastqcDataList.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
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
#' # Return a data_frame/tibble with the raw information
#' Basic_Statistics(fdl)
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Basic_Statistics
#' @aliases Basic_Statistics
setMethod("Basic_Statistics", "FastqcData", function(object){object@Basic_Statistics})

#' @export
#' @rdname Basic_Statistics
#' @aliases Basic_Statistics
setMethod("Basic_Statistics", "FastqcDataList",
          function(object){
            x <- lapply(object@.Data, Basic_Statistics)
            dplyr::bind_rows(x)
          })

#' @export
#' @rdname Basic_Statistics
#' @aliases Basic_Statistics
setMethod("Basic_Statistics", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Basic_Statistics(object)
          })

#' @export
#' @rdname Basic_Statistics
#' @aliases Basic_Statistics
setMethod("Basic_Statistics", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Basic_Statistics(object)
          })

#' @export
#' @rdname Basic_Statistics
#' @aliases Basic_Statistics
setMethod("Basic_Statistics", "character",
          function(object){
            object <- getFastqcData(object)
            Basic_Statistics(object)
          })
