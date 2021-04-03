#' @title Check to see if a file is compressed
#'
#' @description Check to see if a file, or vector of files is compressed
#'
#' @details Reads the first four bytes from the local file header.
#' If the file is a .ZIP file, this should match the magic number
#' `PK\003\004`.
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
#' A `logical` vector
#'
#' @examples
#'
#' # Get the files included with the package
#' fileDir <- system.file("extdata", package = "ngsReports")
#' allFiles <- list.files(fileDir, pattern = "zip$", full.names = TRUE)
#' isCompressed(allFiles)
#'
#' @export
isCompressed <- function(path, type = c("zip", "gzip"), verbose = FALSE){

    stopifnot(file.exists(path))

    ## Define the compression type
    type <- match.arg(type)

    ## This is the magic number for each compression type
    magicNum <- list(zip = c(80, 75, 3, 4), gzip = c(31, 139, 8))[[type]]
    n <- length(magicNum)
    if (verbose > 0) message(
        sprintf("Checking %i file(s) for %s compression", length(path), type)
    )

    ## If any elements of path are found to be directories these should be
    ## identified via a message and not tested, but return FALSE.
    isDir <- file.info(path)$isdir
    if (verbose > 0 & any(isDir)) message(
        sum(isDir), " supplied path(s) represent directories not files:\n",
        paste(path[isDir], collapse = "\n")
    )

    ## Test files which were not found to be directories
    chk <- logical(length(path)) # Defaults to FALSE ensuring dirs are FALSE
    chk[!isDir] <- vapply(path[!isDir], function(x){
        ## Get the magic number from the files
        rw <- readBin(x, what = "raw", n = n)
        if (verbose > 1) message(rw) # Useful for error checking
        sum(rw == magicNum) == n
    }, logical(1))

    ## Return the results, noting again that directories will return FALSE
    chk
}
