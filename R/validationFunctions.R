#' @importFrom checkmate testDirectoryExists
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
    # List the files
    subFiles <- basename(unzip(object@path, list = TRUE)$Name)
  }
  else{
    # Check the file is a directory
    chk <- checkmate::testDirectoryExists(object@path)
    if (!chk){
      warning("The supplied file is not a directory")
      return(FALSE)
    }
    subFiles <- list.files(object@path)
  }
  if (any(!c("fastqc_data.txt", "summary.txt") %in% subFiles)) {
    warning("The required files are missing from the supplied path/file")
    return(FALSE)
  }
  # If it checks out up to here, we're good to go
  TRUE
}

validFastqcFileList <- function(object){
  cls <- vapply(object, class, character(1))
  if (!all(cls == "FastqcFile")) return(FALSE)
}

validFastqcDataList <- function(object){
  # This is very rudimentary & may need more thought
  cls <- vapply(object, class, character(1))
  if (!all(cls == "FastqcData")) return(FALSE)
}
