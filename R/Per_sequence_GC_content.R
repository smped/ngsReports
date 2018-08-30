#' @title Get the Per Sequence GC Content information
#'
#' @description Retrieve the Per Sequence GC Content module from one or more FastQC reports
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
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the Per_sequence_GC_content
#' Per_sequence_GC_content(fdl)   
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "FastqcData",
          function(object){
            df <- object@Per_sequence_GC_content
            if(length(df)){
            df$Filename<- fileName(object)
            dplyr::select(df, "Filename", dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_sequence_GC_content)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_GC_content(object)
          })

#' @export
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_GC_content(object)
          })

#' @export
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "character",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_GC_content(object)
          })
