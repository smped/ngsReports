#' @title Get the Per Base Sequence Quality information
#'
#' @description Retrieve the Per Base Sequence Quality module from one or more FastQC reports
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
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr everything
#' @importFrom dplyr select
#'
#' @export
#' @rdname Per_base_sequence_quality
setMethod("Per_base_sequence_quality", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Per_base_sequence_quality,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname Per_base_sequence_quality
setMethod("Per_base_sequence_quality", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_base_sequence_quality)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname Per_base_sequence_quality
setMethod("Per_base_sequence_quality", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Per_base_sequence_quality(object)
          })

#' @export
#' @rdname Per_base_sequence_quality
setMethod("Per_base_sequence_quality", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Per_base_sequence_quality(object)
          })

#' @export
#' @rdname Per_base_sequence_quality
setMethod("Per_base_sequence_quality", "character",
          function(object){
            object <- getFastqcData(object)
            Per_base_sequence_quality(object)
          })
