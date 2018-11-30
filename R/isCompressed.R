#' @title Check to see if a file is compressed
#'
#' @description Check to see if a file, or vector of files is compressed
#'
#' @details Reads the first four bytes from the local file header.
#' If the file is a .ZIP file, this should match the magic number 
#' \code{PK\003\004}.
#'
#' This function assumes that the first thing in a zip archive is the
#' .ZIP entry with the local file header signature.
#' ZIP files containing a self-extracting archive may not exhibit this structure
#' and will return FALSE
#'
#' @param path The path to one or more files
#' @param type The type of compression to check for.
#' Currently only ZIP/GZIP files have been implemented.
#' @param verbose logical/integer Determine the level of output to show as 
#' messages
#'
#' @return
#' A \code{logical} vector
#'
#' @examples
#'
#' # Get the files included with the package
#' fileDir <- system.file("extdata", package = "ngsReports")
#' allFiles <- list.files(fileDir, pattern = "zip$", full.names = TRUE)
#' isCompressed(allFiles)
#'
#' @export
isCompressed <- function(path, type = "zip", verbose = FALSE){
    
    stopifnot(file.exists(path))
    type <- match.arg(type, c("zip", "gzip"))
    if (verbose > 0) message(sprintf("Checking %i file(s) for %s compression",
                                     length(path), type))
    # This is the length of the magic number for each compression type
    n <- c(zip = 4L, gzip = 3L)[type]
    magicNum <- list(zip = c(80, 75, 3, 4),
                     gzip = c(31, 139, 8))[[type]]
    
    vapply(path, function(x){
        # Get the magic number from the files
        rw <- readBin(x, what = "raw", n = n)
        if (verbose > 1) message(rw)
        sum(rw == magicNum) == n
    }, logical(1))
    
}
