#' Get the set of Basic Statistics
#'
#' @description Retrieve the Basic Statistics module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#' @docType methods
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#' @rdname Basic_Statistics
setMethod("Basic_Statistics", "FastqcData", function(object){object@Basic_Statistics})

#' @export
#' @rdname Basic_Statistics
setMethod("Basic_Statistics", "FastqcDataList",
          function(object){
            x <- lapply(object@.Data, Basic_Statistics)
            dplyr::bind_rows(x)
          })

#' @export
#' @rdname Basic_Statistics
setMethod("Basic_Statistics", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Basic_Statistics(object)
          })

#' @export
#' @rdname Basic_Statistics
setMethod("Basic_Statistics", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Basic_Statistics(object)
          })

#' @export
#' @rdname Basic_Statistics
setMethod("Basic_Statistics", "character",
          function(object){
            object <- getFastqcData(object)
            Basic_Statistics(object)
          })
