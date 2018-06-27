#' @title Get the Sequence Length Distribution information
#'
#' @description Retrieve the Sequence Length Distribution module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#'@examples 
#' 
#' # Get the files included with the package
#' fileDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(fileDir, pattern = "zip$", full.names = TRUE)
#' 
#' # Form a FastqcDataList
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the Sequence Length Distribution
#' Sequence_Length_Distribution(fdl)
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcData",
          function(object){
            df <- object@Sequence_Length_Distribution
            if(length(df)){
            df$Filename <- fileName(object)
            dplyr::select(df, Filename, dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Sequence_Length_Distribution)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Sequence_Length_Distribution(object)
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Sequence_Length_Distribution(object)
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "character",
          function(object){
            object <- getFastqcData(object)
            Sequence_Length_Distribution(object)
          })
