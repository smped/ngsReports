#' @importFrom checkmate testDirectoryExists
isValidFastqcFile <- function(object){

  path <- path(object)
  if (length(path) != 1) {
    warning("Only a single path can be specified")
    return(FALSE)
  }

  comp <- isCompressed(path, type = "zip") # This includes a file.exists step
  if (comp){
    # List the files
    subFiles <- basename(unzip(path, list = TRUE)$Name)
  }
  else{
    # Check the file is a directory
    chk <- checkmate::testDirectoryExists(path)
    if (!chk){
      warning("The supplied file is not a directory")
      return(FALSE)
    }
    subFiles <- list.files(path)
  }
  if (any(!c("fastqc_data.txt", "summary.txt") %in% subFiles)) {
    warning("The required files are missing from the supplied path/file")
    return(FALSE)
  }
  # If it checks out up to here, we're good to go
  TRUE
}

isValidFastqcFileList <- function(object){
  cls <- vapply(object, class, character(1))
  if (!all(cls == "FastqcFile")) return(FALSE)
}

isValidFastqcDataList <- function(object){
  # This is very rudimentary & may need more thought
  cls <- vapply(object, class, character(1))
  if (!all(cls == "FastqcData")) return(FALSE)
}

isValidPwf <- function(object){

  vals <- getColours(object)
  if (length(vals) != 4) return(FALSE) # Ensure all have been set
  if (any(substr(vals, 1, 1) != "#")) return(FALSE) # Start with #
  if (any(!nchar(vals) %in% c(7, 9))) return(FALSE) # RGB length 7 or 9
  if (any(grepl("[G-Z]", vals))) return(FALSE) # Invald Hex values

  TRUE

}
