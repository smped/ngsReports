#' @title Read Hisat2/Bowtie2 Log Files
#'
#' @description Import one or more hisat2/bowtie2 log files as a data frame.
#'
#' @param x \code{character}. Vector of paths to log files
#' @param useBasename \code{logical}. Strip file paths from the Filename column
#' in the returned tibble?
#'
#' @return A \code{tibble}
#'
#' @examples
#' hisat2Logs <- system.file("extdata", "hisat2PE.log", package = "ngsReports")
#' df <- importHisat2Logs(hisat2Logs)
#'
#' # bowtie2 logs are the identical format
#' fl <- c("bowtie2PE.log", "bowtie2SE.log")
#' bowtie2Logs <- system.file("extdata", fl, package = "ngsReports")
#' df <- importBowtie2Logs(bowtie2Logs)
#'
#' @export
importHisat2Logs <- function(x, useBasename = TRUE){

    ## Exit if all files don't exist & let the user get their paths right
    fExists <- file.exists(x)
    if (!all(fExists)) stop(paste("Couldn't find", x[!fExists]))

    data <- lapply(x, readLines)
    validLogs <- vapply(data, .isValidHisat2Log, logical(1))
    if (any(!validLogs)) stop("Incorrect File structure for:\n", x[!validLogs])

    .parseHisat2Data <- function(x){
        x <- stringr::str_trim(x)
        paired <- grepl("were paired", x[2])
        totReads <- gsub("([0-9]*) reads; of these:", "\\1", x[1])
        ## Now find each category then extract the numbers
        noAln <- x[grep("aligned 0 times$", x)]
        noAln <- gsub("([0-9]*) .+ aligned 0 times", "\\1", noAln)
        uniq <- x[grep("aligned exactly 1 time$", x)]
        uniq <- gsub("([0-9]*) .+ aligned exactly 1 time", "\\1", uniq)
        mult <- x[grep("aligned >1 times", x)]
        mult <- gsub("([0-9]*) .+ aligned >1 times", "\\1", mult)

        ## Get the paired only fields
        pairReads <- uniqPairs <- multPairs <- uniqDiscord <- 0
        if (paired) {
            pairReads <- x[grep("were paired; of these:$", x)]
            pairReads <- gsub("([0-9]*) .+ of these:", "\\1", pairReads)
            uniqPairs <- x[grep("aligned concordantly exactly 1 time$", x)]
            uniqPairs <- gsub("([0-9]*) .+ exactly 1 time", "\\1", uniqPairs)
            multPairs <- x[grep("aligned concordantly >1 times$", x)]
            multPairs <- gsub("([0-9]*) .+ >1 times", "\\1", multPairs)
            uniqDiscord <- x[grep("aligned discordantly 1 time$", x)]
            uniqDiscord <- gsub("([0-9]*) .+ 1 time", "\\1", uniqDiscord)
        }

        df <- list(
            Total_Reads = totReads,
            Not_Aligned = noAln,
            Unique_Unpaired = uniq,
            Multiple_Unpaired = mult,
            Paired_Reads = pairReads,
            Unique_In_Pairs = uniqPairs,
            Multiple_In_Pairs = multPairs,
            Unique_Discordant_Pairs = uniqDiscord
        )
        df <- lapply(df, as.integer)
        as_tibble(df)

    }

    out <- lapply(data, .parseHisat2Data)
    out <- dplyr::bind_rows(out)
    out$Filename <- x
    if (useBasename) out$Filename <- basename(out$Filename)
    out$Alignment_Rate <-
        1 - out$Not_Aligned / (out$Total_Reads + out$Paired_Reads)

    ## Now return in the correct order
    dplyr::select(out,
                  "Filename",
                  tidyselect::ends_with("Reads"),
                  tidyselect::contains("Unique"),
                  tidyselect::contains("Multiple"),
                  tidyselect::everything())

}

#' @rdname importHisat2Logs
#' @aliases importBowtie2Logs
#' @export
importBowtie2Logs <- importHisat2Logs

#' @title Check for a valid Hisat2 log
#' @description Checks internal structure of the parsed log file
#' @param x Character vector containing parsed log file using the function
#' readLines
#' @return logical(1)
#' @keywords internal
.isValidHisat2Log <- function(x){
    n <- length(x)
    chkLen <- length(x) > 0
    firstLine <- grepl("reads; of these:$", x[1])
    lastLine <- grepl("overall alignment rate$", x[n])
    noAln <- sum(grepl("aligned 0 times$", x)) == 1
    alnExact <- sum(grepl("aligned exactly 1 time$", x)) == 1
    alnG1 <- sum(grepl("aligned >1 times$", x)) == 1
    all(c(chkLen, firstLine, lastLine, noAln, alnExact, alnG1))
}
