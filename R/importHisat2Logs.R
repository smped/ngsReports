#' @title Read Hisat2/Bowtie2 Log Files
#'
#' @description Import one or more hisat2/bowtie2 log files as a data frame.
#'
#' @param x \code{character}. Vector of paths to log files
#'
#' @return A \code{data_frame}
#'
#' @examples
#' hisat2Logs <- system.file("extdata", "hisat2PE.log", package = "ngsReports")
#' df <- importHisat2Logs(hisat2Logs)
#'
#' # bowtie2 logs are the identical format
#' bowtie2Logs <- system.file("extdata", c("bowtie2PE.log", "bowtie2SE.log"), package = "ngsReports")
#' df <- importBowtie2Logs(bowtie2Logs)
#'
#' @export
importHisat2Logs <- function(x){

  fe <- file.exists(x)
  if (all(!fe)) stop("Files could not be found")
  if (any(!fe) ){
    warning("The following file(s) could not be found:\n", paste(x[!fe], sep = "\n"),
            "\nThese files will be ignored")
  }
  x <- x[fe]

  data <- lapply(x, readLines)
  names(data) <- basename(x)

  isHisat2Log <- function(x){
    n <- length(x)
    if (n == 0) return(FALSE)
    if (!grepl("reads; of these:", x[1])) return(FALSE)
    if (!grepl("overall alignment rate", x[n])) return(FALSE)
    TRUE
  }

  getHisat2Data <- function(x){
    paired <- grepl("were paired", x[2])
    df <- data.frame(Total_Reads = as.integer(gsub("([0-9]*) reads; of these:", "\\1", x[1])))
    df$Not_Aligned <- as.integer(gsub("([0-9]*) \\(.+\\) aligned 0 times", "\\1",
                                      x[grep("aligned 0 times$", x)]))
    df$Unique_Unpaired <- as.integer(gsub("([0-9]*) \\(.+\\) aligned exactly 1 time", "\\1",
                                          x[grep("aligned exactly 1 time", x)]))
    df$Multiple_Unpaired <- as.integer(gsub("([0-9]*) \\(.+\\) aligned >1 times", "\\1",
                                            x[grep("aligned >1 times", x)]))
    if (paired){
      df$Paired_Reads <- as.integer(gsub("([0-9]*) \\(.+\\) were paired; of these:", "\\1",
                                         x[grep("were paired; of these:", x)]))
      df$Unique_In_Pairs <- as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly exactly 1 time", "\\1",
                                            x[grep("aligned concordantly exactly 1 time", x)]))
      df$Multiple_In_Pairs <- as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly >1 times", "\\1",
                                              x[grep("aligned concordantly >1 times", x)]))
      df$Unique_Discordant_Pairs <- as.integer(gsub("([0-9]*) \\(.+\\) aligned discordantly 1 time", "\\1",
                                                    x[grep("aligned discordantly 1 time", x)]))
      df$Alignment_Rate = with(df, 1 - Not_Aligned / (Total_Reads + Paired_Reads))
    }
    else{
      df$Alignment_Rate = with(df, 1 - Not_Aligned / Total_Reads)
    }
    df
  }

  validLogs <- vapply(data, isHisat2Log, logical(1))
  if (any(!validLogs)) stop("Incorrect File structure for:\n", names(validLogs)[!validLogs])

  out <- lapply(data, getHisat2Data)
  out <- dplyr::bind_rows(out)
  out$Filename <- names(data)

  dplyr::select(out, Filename,
                dplyr::ends_with("Reads"),
                dplyr::contains("Unique"),
                dplyr::contains("Multiple"),
                dplyr::everything())

}

#' @rdname importHisat2Logs
#' @aliases importBowtie2Logs
#' @export
importBowtie2Logs <- importHisat2Logs
