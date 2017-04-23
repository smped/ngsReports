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
#' @include AllGenerics.R
#'
#' @export
#' @rdname FastqcFile-methods
#' @aliases FastqcFile,character-method
setMethod("FastqcFile", "character",
          function(filePath){
            # Read the first 4 bytes as hexadecimal values
            rw <- readBin(filePath, what = "raw", n = 4L)
            # Zipped files start with c(80, 75, 03, 04) in the first 4 bytes
            zipBytes <- as.raw(c(80, 75, 03, 04))
            # Comparison below will give 4 TRUE values for a zipped file
            # If filePath is a directory, it will return raw(0)
            # which gives logical(0) when compared to `zipBytes`, hence the sum == 4
            comp <- sum(rw == zipBytes) == 4
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
