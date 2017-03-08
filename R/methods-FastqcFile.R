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
#' @include AllClasses.R
#'
#' @export
#' @rdname FastqcFile-methods
#' @aliases FastqcFile,character-method
setMethod("FastqcFile", "character",
          function(filePath){
            comp <- grepl("zip$", filePath)
            new("FastqcFile", path = filePath, compressed = comp)
          })

#' @export
#' @rdname FastqcFile-methods
#' @aliases path,FastqcFile-method
setMethod("path", "FastqcFile", function(object){object@path})

#' @export
#' @rdname FastqcFile-methods
#' @aliases isCompressed,FastqcFile-method
setMethod("isCompressed", "FastqcFile", function(object){object@compressed})

#' @export
#' @rdname FastqcFile-methods
#' @aliases names,FastqcFile-method
setMethod("names", "FastqcFile", function(x){basename(x@path)})

#' @export
#' @rdname FastqcFile-methods
#' @aliases show,FastqcFile-methods
setMethod("show", "FastqcFile",
          function(object){
            cat(names(object), "\n")
            cat("Located in", dirname(path(object)), "\n")
          })
