#' @title Get the Overrepresented Sequences information
#'
#' @description Retrieve the Overrepresented Sequences module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#'
#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Overrepresented_sequences)
            if(length(df)){
            df$Filename <- fileName(object)
            dplyr::select(df, Filename, dplyr::everything())
            }
            else NULL
          })

#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Overrepresented_sequences)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Overrepresented_sequences(object)
          })

#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Overrepresented_sequences(object)
          })

#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "character",
          function(object){
            object <- getFastqcData(object)
            Overrepresented_sequences(object)
          })
