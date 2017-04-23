#' @title Get the Sequence Length Distribution information
#'
#' @description Retrieve the Sequence Length Distribution module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution,FastqcData-method
setMethod("Sequence_Length_Distribution", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Sequence_Length_Distribution,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })
#'
#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution,FastqcDataList-method
setMethod("Sequence_Length_Distribution", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Sequence_Length_Distribution)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Sequence_Length_Distribution(object)
          })

#' @export
#' @rdname Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Sequence_Length_Distribution(object)
          })

#' @export
#' @rdname Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "character",
          function(object){
            object <- getFastqcData(object)
            Sequence_Length_Distribution(object)
          })
