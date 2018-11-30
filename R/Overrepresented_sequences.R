#' @title Get the Overrepresented Sequences information
#'
#' @description Retrieve the Overrepresented Sequences module from one or more 
#' FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, 
#' \code{FastqcData}, \code{fastqcDataList}, or simply a \code{character} vector
#' of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{tibble} containing all information combined from all 
#' supplied FastQC reports
#' @examples 
#' 
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the Overrepresented Sequences
#' Overrepresented_sequences(fdl)
#'
#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "FastqcData",
          function(object){
              df <- object@Overrepresented_sequences
              if (length(df)) {# Check there is data in the module
                  # Add a Filename column if there is any data
                  df$Filename <- fileName(object)
                  dplyr::select(df, "Filename", tidyselect::everything())
              }
              else {# Otherwise return the blank data.frame
                  df
              }
          })

#' @export
#' @rdname Overrepresented_sequences
#' @aliases Overrepresented_sequences
setMethod("Overrepresented_sequences", "FastqcDataList",
          function(object){
              df <- lapply(object@.Data, Overrepresented_sequences)
              nulls <- vapply(df, 
                              function(x){
                                  length(x) == 0
                              }, logical(1))
              if (sum(nulls) > 0) message(
                  sprintf("The Overrepresented_sequences module was empty in:\n%s",
                          paste(path(object)[nulls], sep = "\n"))
              )
              dplyr::bind_rows(df)
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
