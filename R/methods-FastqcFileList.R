#' Path to FastQC data a vector of fastq files
#'
#' @description Define a FastqcFileList
#'
#' @details Checks the structure of a folder or zip file with the output from FastQC
#'
#' @param path The path to a vector of zip files (or uncompressed folders) as output by FastQC
#' @param object An object of class \code{FastqcFileList}
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' 
#'
#' @export
setMethod("FastqcFileList", "character",
          function(path)
          {
            fls <- lapply(path, FastqcFile)
            names(fls) <- basename(path)
            new("FastqcFileList", fls)
          })

#' @export
setMethod("FastqcFileList", "list",
          function(path)
          {
            cls <- vapply(path, class, character(1))
            if (any(!cls %in% "FastqcFile")) stop("Method can only be applied to\nFastqcFile objects as a generic list.")
            new("FastqcFileList", path)
          })

#' @export
setMethod("path", "FastqcFileList",
          function(object){vapply(object, path, character(1))})

#' @export
setMethod("fileNames", "FastqcFileList", function(object){vapply(object, fileNames, character(1))})

# Define the subsetting method, to remain consistent with list behaviour
setMethod("[", "FastqcFileList", function(x, i, j, ..., drop = TRUE){FastqcFileList(x@.Data[i])})

# The show method doesn't need exporting
setMethod("show", "FastqcFileList",
          function(object){
            l <- length(object)
            # cmp <- sum(isCompressed(object))
            cat("FastqcFileList of", l, "file(s).\n")
            cat("Located in:\n", paste(unique(dirname(path(object))), collapse = "\n"))
          })
