#' @title Get the Kmer Content information
#'
#' @description Retrieve the Kmer Content module from one or more FastQC reports
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
#' # Print the Kmer_Content
#' Kmer_Content(fdl)
#'
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcData",
          function(object){
            df <- dplyr::mutate(object@Kmer_Content)
            if(length(df)){
            df$Filename <- fileName(object)
            dplyr::select(df, "Filename", dplyr::everything())
            }
            else NULL
            })

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, Kmer_Content)
            if(length(unlist(df))) dplyr::bind_rows(df)
            else NULL
          })

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcFile",
          function(object){
            object <- getFastqcData(object)
            Kmer_Content(object)
          })

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcFileList",
          function(object){
            object <- getFastqcData(object)
            Kmer_Content(object)
          })

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "character",
          function(object){
            object <- getFastqcData(object)
            Kmer_Content(object)
          })
