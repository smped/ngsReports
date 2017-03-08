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

validFastqcFileList <- function(object){
  cls <- vapply(object, class, character(1))
  all(cls == "FastqcFile")
}
