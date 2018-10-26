#' @title Get the Sequence Length Distribution information
#'
#' @description Retrieve the Sequence Length Distribution module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{tibble} containing all information combined from all supplied FastQC reports
#'
#'@examples 
#' 
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the Sequence Length Distribution
#' Sequence_Length_Distribution(fdl)
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcData",
          function(object){
              df <- object@Sequence_Length_Distribution
              if(length(df)){ # Check there is data in the module
                  # Add a Filename column if there is any data
                  df$Filename <- fileName(object)
                  dplyr::select(df, "Filename", tidyselect::everything())
              }
              else { # Otherwise return the blank data.frame
                  df
              }
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcDataList",
          function(object){
              df <- lapply(object@.Data, Sequence_Length_Distribution)
              nulls <- vapply(df, 
                              function(x){
                                  length(x) == 0
                              }, logical(1))
              if (sum(nulls) > 0) message(
                  sprintf("The Sequence_Length_Distribution module was missing from:\n%s",
                          paste(path(object)[nulls], sep = "\n"))
              )
              dplyr::bind_rows(df)
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcFile",
          function(object){
              object <- getFastqcData(object)
              Sequence_Length_Distribution(object)
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "FastqcFileList",
          function(object){
              object <- getFastqcData(object)
              Sequence_Length_Distribution(object)
          })

#' @export
#' @rdname Sequence_Length_Distribution
#' @aliases Sequence_Length_Distribution
setMethod("Sequence_Length_Distribution", "character",
          function(object){
              object <- getFastqcData(object)
              Sequence_Length_Distribution(object)
          })
