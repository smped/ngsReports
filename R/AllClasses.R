#' Define object classes
#'
#' @slot path a character representation of a path to a single file/directory
#' @slot compressed logical. Is the FastQC file compressed as a zip file
FastqcFile <- setClass("FastqcFile",
                       slots = c(path = "character",
                                 compressed = "logical"))

FastqcFileList <- setClass("FastqcFileList", contains="FastqcFile")

validFastqcFile <- function(object){

  if (length(object@path) != 1) {
    warning("Only a single path can be specified")
    return(FALSE)
  }

  if (!file.exists(object@path)) {
    warning("Supplied path not found")
    return(FALSE)
  }

  if (object@compressed){
    # Check the file is actually compressed
    if(!grepl("zip$", object@path)) {
      warning("Path does not point to a compressed file")
      return(FALSE)
    }
    # List the files
    subFiles <- basename(unzip(object@path, list = TRUE)$Name)
  }
  else{
    subFiles <- list.files(object@path)
  }
  if (any(!c("fastqc_data.txt", "summary.txt") %in% subFiles)) {
    warning("The required files are missing from the supplied path")
    return(FALSE)
  }
  # If it checks out up to here, we're good to go
  TRUE
}

setValidity("FastqcFile", validFastqcFile)
