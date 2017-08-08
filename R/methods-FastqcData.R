#' Methods for an object of class FastqcData
#'
#' @description Methods and Accessors for an object of class FastqcData
#'
#' @details Describes all methods for an object of class \code{FastqcData}
#'
#' @param object An object of class \code{FastqcData}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' 
#'
#' @export
setMethod("path", "FastqcData", function(object){object@path})

#' @export
setMethod("Version", "FastqcData", function(object){object@Version})

#' @export
setMethod("getSummary", "FastqcData", function(object){object@Summary})

#' @export
setMethod("fileNames", "FastqcData", function(object){object@Summary$Filename[1]})

# The show method doesn't need exporting
setMethod("show", "FastqcData",
          function(object){
            cat("FastqcData for", object@Basic_Statistics$Filename, "\n")
            cat("Source Fastq file contains", scales::comma(object@Basic_Statistics$Total_Sequences), "reads.\n")
            cat("Source FastQC file is located in", object@path)
          })
