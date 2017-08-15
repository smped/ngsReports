#' @title Read Hisat2 Log Files
#'
#' @description Import one or more hisat2 log files as a data frame.
#'
#' @param logs \code{character}. Vector of paths to hisat2 log files
#'
#' @return A \code{data_frame}
#'
#' @export
importHisat2Logs <- function(logs){

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
    n <- length(x)
    if (n == 0) return(FALSE)
    if (!grepl("reads; of these:", x[1])) return(FALSE)
    if (!grepl("overall alignment rate", x[n])) return(FALSE)
    TRUE
  }

  getHisat2Data <- function(x){
    paired <- grepl("were paired", x[2])
    if (paired){
      df <- dplyr::data_frame(
        Total_Reads = as.integer(gsub("([0-9]*) reads; of these:", "\\1", x[1])),
        Paired_Reads = as.integer(gsub("([0-9]*) \\(.+\\) were paired; of these:", "\\1", x[2])),
        Unique_In_Pairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly exactly 1 time", "\\1", x[4])),
        Multiple_In_Pairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly >1 times", "\\1", x[5])),
        Unique_Discordant_Pairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned discordantly 1 time", "\\1", x[8])),
        Unique_Unpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned exactly 1 time", "\\1", x[13])),
        Multiple_Unpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned >1 times", "\\1", x[14])),
        Not_Aligned = as.integer(gsub("([0-9]*) \\(.+\\) aligned 0 times", "\\1", x[12])))
      df$Alignment_Rate = with(df, 1 - Not_Aligned / (Total_Reads + Paired_Reads))
    }
    else{
      df <- dplyr::data_frame(
        Total_Reads = as.integer(gsub("([0-9]*) reads; of these:", "\\1", x[1])),
        Not_Aligned = as.integer(gsub("([0-9]*) \\(.+\\) aligned 0 times", "\\1", x[3])),
        Unique_Unpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned exactly 1 time", "\\1", x[4])),
        Multiple_Unpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned >1 times", "\\1", x[5])))
      df$Alignment_Rate = with(df, 1 - Not_Aligned / Total_Reads)
    }
    df
  }

  validLogs <- vapply(data, isHisat2Log, logical(1))
  if (any(!validLogs)) warning("Incorrect File structure for:\n", names(validLogs)[!validLogs])
  data <- data[validLogs]

  out <- lapply(data, getHisat2Data)
  out <- dplyr::bind_rows(out)
  out$Filename <- names(data)

  dplyr::select(out, Filename, dplyr::everything())

}
