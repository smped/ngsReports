#' @title Read Hisat2 Log Files
#'
#' @description Import one or more hisat2 log files as a data frame.
#'
#' @param logs \code{character}. Vector of paths to hisat2 log files
#'
#' @return A \code{data_frame}
#'
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#' @importFrom dplyr everything
#'
#' @export
readHisat2Logs <- function(logs){

  fe <- file.exists(logs)
  if (all(!fe)) stop("Files could not be found")
  if (any(!fe) ){
    warning("The following file(s) could not be found:\n", paste(logs[!fe], sep = "\n"),
            "\nThese files will be ignored")
  }
  logs <- logs[fe]

  data <- lapply(logs, readLines)
  names(data) <- basename(logs)

  isHisat2Log <- function(x){
    if (length(x) != 15) return(FALSE)
    if (grep("reads; of these:", x) != 1) return(FALSE)
    if (any(as.logical(grep("concordantly", x) - c(3:5, 7, 10)))) return(FALSE)
    if (any(as.logical(grep("discordantly", x) - c(8, 10)))) return(FALSE)
    if (grep("overall alignment rate", x) != 15) return(FALSE)
    TRUE
  }

  getHisat2Data <- function(x){
    dplyr::data_frame(
      TotalReads = as.integer(gsub("([0-9]*) reads; of these:", "\\1", x[1])),
      PairedReads = as.integer(gsub("([0-9]*) \\(.+\\) were paired; of these:", "\\1", x[2])),
      UniqueInPairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly exactly 1 time", "\\1", x[4])),
      MultipleInPairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly >1 times", "\\1", x[5])),
      UniqueDiscordantPairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned discordantly 1 time", "\\1", x[8])),
      UniqueUnpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned exactly 1 time", "\\1", x[13])),
      MultipleUnpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned >1 times", "\\1", x[14])),
      NotAligned = as.integer(gsub("([0-9]*) \\(.+\\) aligned 0 times", "\\1", x[12])),
      AlignmentRate = 1 - NotAligned / (TotalReads + PairedReads)
    )
  }

  validLogs <- vapply(data, isHisat2Log, logical(1))
  if (any(!validLogs)) warning("Incorrect File structure for:\n", names(validLogs)[!validLogs])
  data <- data[validLogs]

  out <- lapply(data, getHisat2Data)
  out <- dplyr::bind_rows(out)
  out$Filename <- names(data)

  dplyr::select(out, Filename, dplyr::everything())

}
