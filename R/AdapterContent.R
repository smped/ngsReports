#' Get the Adapter Content information
#'
#' @description Retrieve the Adapter Content module from one or more FastQC reports
#'
#' @details It is currently assumed that all FastQC reports contain information about the same adapters.
#' If this was varied during the preparation of the reports, this function will probably fail.
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files.
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
#' Adapter_Content(fdl)
#'
#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod("Adapter_Content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Adapter_Content)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod("Adapter_Content", "FastqcData",
          function(object){
            df <- object@Adapter_Content
            if(length(df)){
              df$Filename <- fileName(object)
              dplyr::select(df, "Filename", dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod("Adapter_Content", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Adapter_Content(object)
          })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod("Adapter_Content", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Adapter_Content(object)
          })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod("Adapter_Content", "character",
          function(object){
            object <- getFastqcData(object)
            Adapter_Content(object)
          })
