#' Load STAR Alignment Summaries
#'
#' @description Loads the information from one or more Log.final.out files as output by STAR
#'
#' @param path \code{character}. Path to one or more STAR alignment summaries
#'
#' @return Returns an object of class \code{starLogData}
#'
#' @importFrom readr read_delim
#' @importFrom lubridate parse_date_time
#' @importFrom dplyr data_frame
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#'
#' @export
readStarLogs <- function(path){

  path <- unique(path) # Remove any  duplicates
  stopifnot(file.exists(path)) # Check they all exist
  # nFiles <- length(path)
  data <- suppressWarnings(lapply(path, readr::read_delim,
                 delim = "\t", col_names = c("Category", "Value"),
                 col_types = "cc", trim_ws = TRUE) )# Load them all
  names(data) <- basename(path)

  isValidStarLog <- function(x){
    nLines <- length(x)
    if (!grepl("Started job on", x$Category[1])) return(FALSE)
    if (!grepl("UNIQUE READS:", x$Category[7])) return(FALSE)
    if (!grepl("MULTI-MAPPING READS:", x$Category[22])) return(FALSE)
    if (!grepl("UNMAPPED READS:", x$Category[27])) return(FALSE)
    TRUE
  }

  validLogs <- vapply(data, isValidStarLog, logical(1))
  if (any(!validLogs)) {
    stop(paste("Incorrect file structure for:", names(validLogs)[!validLogs], collapse = "\n"))
  }

  getTiming <- function(x){
    times <- lubridate::parse_date_time(x$Value[1:3], orders = "b! d! HMS")
    dplyr::data_frame(
      StartedJob = times[1],
      StartedMapping = times[2],
      Finished = times[3],
      TotalTime = difftime(Finished, StartedJob),
      MappingTime = difftime(Finished, StartedMapping),
      Speed = filter(x, grepl("Mapping speed", Category))$Value
      )
  }

  getInput <- function(x){
    vals <- as.list(dplyr::filter(x, grepl("input", Category))$Value)
    names(vals) <- c("TotalReads", "AverageLength")
    vals <- lapply(vals, as.integer)
    dplyr::as_data_frame(vals)
  }

  getUnique <- function(x){
    dplyr::data_frame(
      UniquelyMapped = as.integer(filter(x, grepl("Uniquely mapped reads number", Category))$Value),
      UniqueMappingRate = UniquelyMapped / getInput(x)$TotalReads,
      AveMappedLength = as.double(filter(x, grepl("Average mapped length", Category))$Value)
    )
  }

  getSplicing <- function(x){
    vals <- dplyr::filter(x, grepl("splices", Category))$Value
    vals <- as.integer(vals)
    names(vals) <- c("Total", "Annotated", "GT/AG", "GC/AG", "AT/AC", "Non-canonical")
    dplyr::as_data_frame(as.list(vals))
  }

  getMulti <- function(x){
    vals <- dplyr::filter(x, grepl("Number of reads.+loci", Category))$Value
    vals <- as.integer(vals)
    names(vals) <- c("MultipleLoci", "TooManyLoci")
    df <- dplyr::as_data_frame(as.list(vals))
    df <- mutate(df,
                 MultiMappingRate = MultipleLoci / getInput(x)$TotalReads,
                 TooManyRate = TooManyLoci / getInput(x)$TotalReads)
    dplyr::select(df, contains("Multi"), contains("TooMany"))
  }

  getErrors <- function(x){
    vals <- dplyr::filter(x, grepl("(Mismatch|Deletion|Insertion)", Category))$Value
    names(vals) <- c("MismatchRate", "DeletionRate", "AveDeletionLength",
                     "InsertionRate", "AveInsertionLength")
    vals <- as.list(vals)
    vals$MismatchRate = as.numeric(gsub("\\%", "", vals$MismatchRate)) / 100
    vals$DeletionRate = as.numeric(gsub("\\%", "", vals$DeletionRate)) / 100
    vals$InsertionRate = as.numeric(gsub("\\%", "", vals$InsertionRate)) / 100
    dplyr::as_data_frame(vals)
  }

  getUnmapped <- function(x){
    x <- dplyr::filter(x, grepl("unmapped", Category))
    vals <- x$Value
    vals <- as.numeric(gsub("\\%", "", vals)) / 100
    nm <- gsub(".*\\% of reads unmapped: (.+) *\\|.*", "\\1", x$Category)
    nm <- stringr::str_trim(nm)
    nm <- stringr::str_to_title(nm)
    names(vals) <- nm
    dplyr::as_data_frame(as.list(vals))
  }

  new("StarLogData",
      fileName = basename(path),
      timing = lapply(data, getTiming) %>% dplyr::bind_rows(),
      input = lapply(data, getInput) %>% dplyr::bind_rows(),
      unique = lapply(data, getUnique) %>% dplyr::bind_rows(),
      splicing = lapply(data, getSplicing) %>% dplyr::bind_rows(),
      errors  = lapply(data, getErrors) %>% dplyr::bind_rows(),
      multiMapping = lapply(data, getMulti) %>% dplyr::bind_rows(),
      unmapped = lapply(data, getUnmapped) %>% dplyr::bind_rows()
  )


}
