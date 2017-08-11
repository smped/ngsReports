#' @title Get the Per Sequence GC Content information
#'
#' @description Retrieve the Per Sequence GC Content module from one or more FastQC reports
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
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "FastqcData",
          function(object){
            df <- object@Per_sequence_GC_content
            df$Filename<- fileNames(object)
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname Per_sequence_GC_content
#' @aliases Per_sequence_GC_content
setMethod("Per_sequence_GC_content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Per_sequence_GC_content)
            dplyr::bind_rows(df)
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
