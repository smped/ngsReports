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
            # Zipped files start with c(80, 75, 03, 04) in the first 4 bytes
            rw <- readBin(filePath, what = "raw", n = 4L)
            comp <- sum(rw == as.raw(c(80, 75, 03, 04))) == 4
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
#' @aliases fileNames,FastqcFile-method
setMethod("fileNames", "FastqcFile", function(object){basename(object@path)})

#' @export
#' @rdname FastqcFile-methods
#' @aliases show,FastqcFile-methods
setMethod("show", "FastqcFile",
          function(object){
            cat(names(object), "\n")
            cat("Located in", dirname(path(object)), "\n")
          })
