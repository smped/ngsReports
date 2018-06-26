#' @title Get the Sequence Duplication Levels information
#'
#' @description Retrieve the Sequence Duplication Levels module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcData",
          function(object){
            df <- object@Sequence_Duplication_Levels
            if(length(df)){
            df$Filename <- fileName(object)
            dplyr::select(df, Filename, dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Sequence_Duplication_Levels)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Sequence_Duplication_Levels(object)
          })

#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Sequence_Duplication_Levels(object)
          })

#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "character",
          function(object){
            object <- getFastqcData(object)
            Sequence_Duplication_Levels(object)
          })
