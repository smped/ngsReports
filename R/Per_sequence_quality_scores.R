#' @title Get the Per Sequence Quality Scores information
#'
#' @description Retrieve the Per Sequence Quality Scores module from one or more FastQC reports
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
#' @export
#' @rdname Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_sequence_quality_scores,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_sequence_quality_scores)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_quality_scores(object)
          })

#' @export
#' @rdname Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_quality_scores(object)
          })

#' @export
#' @rdname Per_sequence_quality_scores
setMethod("Per_sequence_quality_scores", "character",
          function(object){
            object <- getFastqcData(object)
            Per_sequence_quality_scores(object)
          })
