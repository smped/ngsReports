#' Path to FastQC data from a single fastq file
#'
#' @description Define a FastqcFile
#'
#' @details Checks the structure of a folder or zip file with the output from FastQC
#'
#' @param filePath The path to a single FastQC output zip file, or uncompressed folder
#' @param object An object of class \code{FastqcFile}
#'
#' @return An object of cass \code{FastqcFile}
#'
#' @docType methods
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @export
setMethod("FastqcFile", "character",
          function(filePath){
            new("FastqcFile", path = filePath)
          })

#' @export
setMethod("path", "FastqcFile", function(object){object@path})

#' @export
setMethod("fileNames", "FastqcFile", function(object){basename(object@path)})

# The show method doesn't need exporting
setMethod("show", "FastqcFile",
          function(object){
            cat(fileNames(object), "\n")
            cat("Located in", dirname(path(object)), "\n")
          })
