#' Path to data from FastQC
#'
#' @description Define a FastqcFile
#'
#' @details Checks the structure of a folder or zip file with the output from FastQC
#'
#' @param path The path to a single FastQC output zip file, or uncompressed folder
#'
#' @return An object of cass \code{FastqcFile}
#'
#' @include AllClasses.R
#'
#' @export
#' @docType methods
#' @rdname FastqcFile-methods
setGeneric("FastqcFile",
           function(path)
             standardGeneric("FastqcFile"),
           signature = "path")

#' @rdname FastqcFile-methods
#' @aliases FastqcFile,character-method
setMethod(FastqcFile, "character",
          function(path){
            comp <- grepl("zip$", path)
            new("FastqcFile", path = path, compressed = comp)
          })

#' @export
#' @rdname FastqcFile-methods
#' @aliases FastqcFileList,character-method
setGeneric("FastqcFileList",
           function(path)
             standardGeneric("FastqcFileList"),
           signature="path")

#' @rdname FastqcFile-methods
#' @aliases FastqcFileList,character-method
setMethod(FastqcFileList, "character",
          function(path)
          {
            fls <- lapply(path, FastqcFile)
            class(fls) <- "FastqcFileList"
          })
