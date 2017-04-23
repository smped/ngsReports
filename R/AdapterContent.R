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
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#' @export
#' @rdname Adapter_Content
setMethod("Adapter_Content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Adapter_Content)
            dplyr::bind_rows(df)
          })

#' @export
#' @rdname Adapter_Content
setMethod("Adapter_Content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Adapter_Content,
                                Filename = fileNames(object))
            dplyr::select(df, Filename, dplyr::everything())
          })

#' @export
#' @rdname Adapter_Content
setMethod("Adapter_Content", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Adapter_Content(object)
          })

#' @export
#' @rdname Adapter_Content
setMethod("Adapter_Content", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Adapter_Content(object)
          })

#' @export
#' @rdname Adapter_Content
setMethod("Adapter_Content", "character",
          function(object){
            object <- getFastqcData(object)
            Adapter_Content(object)
          })
