#' @title Get the Per Base Sequence Content information
#'
#' @description Retrieve the Per Base Sequence Content module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData}, \code{fastqcDataList},
#' or simply a \code{character} vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{data_frame} containing all information combined from all supplied FastQC reports
#'
#' @examples 
#' 
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' 
#' # Print the Per_base_sequence_content
#' Per_base_sequence_content(fdl)   
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Per_base_sequence_content
#' @aliases Per_base_sequence_content
setMethod("Per_base_sequence_content", "FastqcData",
          function(object){
              df <- object@Per_base_sequence_content
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
#' @rdname Per_base_sequence_content
#' @aliases Per_base_sequence_content
setMethod("Per_base_sequence_content", "FastqcDataList",
          function(object){
              df <- lapply(object@.Data, Per_base_sequence_content)
              nulls <- vapply(df, 
                              function(x){
                                  length(x) == 0
                              }, logical(1))
              if (sum(nulls) > 0) message(
                  sprintf("The Per_base_sequence_content module was missing from:\n%s",
                          paste(path(object)[nulls], sep = "\n"))
              )
              dplyr::bind_rows(df)
          })

#' @export
#' @rdname Per_base_sequence_content
#' @aliases Per_base_sequence_content
setMethod("Per_base_sequence_content", "FastqcFile",
          function(object){
              object <- getFastqcData(object)
              Per_base_sequence_content(object)
          })

#' @export
#' @rdname Per_base_sequence_content
#' @aliases Per_base_sequence_content
setMethod("Per_base_sequence_content", "FastqcFileList",
          function(object){
              object <- getFastqcData(object)
              Per_base_sequence_content(object)
          })

#' @export
#' @rdname Per_base_sequence_content
#' @aliases Per_base_sequence_content
setMethod("Per_base_sequence_content", "character",
          function(object){
              object <- getFastqcData(object)
              Per_base_sequence_content(object)
          })
