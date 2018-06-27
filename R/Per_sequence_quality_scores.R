#' @title Get the Per Sequence Quality Scores information
#'
#' @description Retrieve the Per Sequence Quality Scores module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#' @examples 
#' 
#' # Get the files included with the package
#' fileDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(fileDir, pattern = "zip$", full.names = TRUE)
#' 
#' # Form a FastqcDataList
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the Per_sequence_quality_scores
#' Per_sequence_quality_scores(fdl)   
#'
#' @docType methods
#'
#' @export
#' @rdname Per_sequence_quality_scores
#' @aliases Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcData",
          function(object){
            df <- object@Per_sequence_quality_scores
            if(length(df)){
            df$Filename<- fileName(object)
            dplyr::select(df, Filename, dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Per_sequence_quality_scores
#' @aliases Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_sequence_quality_scores)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
            })

#' @export
#' @rdname Per_sequence_quality_scores
#' @aliases Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_quality_scores(object)
          })

#' @export
#' @rdname Per_sequence_quality_scores
#' @aliases Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_quality_scores(object)
          })

#' @export
#' @rdname Per_sequence_quality_scores
#' @aliases Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "character",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_quality_scores(object)
          })
