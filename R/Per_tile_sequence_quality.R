#' @title Get the Per Tile Sequence Quality information
#'
#' @description Retrieve the Per Tile Sequence Quality module from one or more FastQC reports
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
#'
#' @export
#' @rdname Per_tile_sequence_quality
#' @aliases Per_tile_sequence_quality
setMethod("Per_tile_sequence_quality", "FastqcData",
          function(object){
            df <- object@Per_tile_sequence_quality
            if(length(df)){
            df$Filename <- fileName(object)
            dplyr::select(df, Filename, dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Per_tile_sequence_quality
#' @aliases Per_tile_sequence_quality
setMethod("Per_tile_sequence_quality", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_tile_sequence_quality)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Per_tile_sequence_quality
#' @aliases Per_tile_sequence_quality
setMethod("Per_tile_sequence_quality", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Per_tile_sequence_quality(object)
          })

#' @export
#' @rdname Per_tile_sequence_quality
#' @aliases Per_tile_sequence_quality
setMethod("Per_tile_sequence_quality", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Per_tile_sequence_quality(object)
          })

#' @export
#' @rdname Per_tile_sequence_quality
#' @aliases Per_tile_sequence_quality
setMethod("Per_tile_sequence_quality", "character",
          function(object){
            object <- getFastqcData(object)
            Per_tile_sequence_quality(object)
          })
