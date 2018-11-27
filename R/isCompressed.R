#' Check to see if a file is compressed
#'
#' @description Check to see if a file, or vector of files is compressed
#'
#' @details Reads the first four bytes from the local file header.
#' If the file is a .ZIP file, this should match the magic number \code{PK\003\004}.
#'
#' This function assumes that the first thing in a zip archive is the
#' .ZIP entry with the local file header signature.
#' ZIP files containing a self-extracting archive may not exhibit this structure
#' and will return FALSE
#'
#' @param path The path to one or more files
#' @param type The type of compression to check for.
#' Currently only ZIP files have been implemented.
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
isCompressed <- function(path, type = "zip"){

  stopifnot(file.exists(path))
  if (type != "zip") stop("Currently only zip files are implemented")

  vapply(path, function(x){
    # Read the first 4 bytes as hexadecimal values
    rw <- readBin(x, what = "raw", n = 4L)
    rawToChar(rw) == "PK\003\004" # Magic number
  }, logical(1))

}
