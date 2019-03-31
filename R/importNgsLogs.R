#' @title Import Various NGS-related log files
#'
#' @description Imports NGS-related log files such as those generated from
#' stderr
#'
#' @details Imports one or more log files as output by tools such as:
#' \code{bowtie}, \code{bowtie2}, \code{Hisat2} or \code{STAR}.
#'
#' @param x \code{character}. Vector of filenames. All log files must be of the
#' same type. Duplicate file paths will be silently ignored.
#' @param type \code{character}. The type of file being imported. Can be one of
#' \code{bowtie}, \code{bowtie2}, \code{hisat2} or \code{star}
#' @param which Which element of the parsed object to return. Ignored in all
#' file types except duplication metrics which can return either the metrics
#' or the data supplied as a histogram. Defaults to the metrics data in this
#' case.
#'
#' @return A \code{tibble}.
#' Column names are broadly similar to the text in supplied files,
#' but have been modified for easier handling under R naming conventions.
#'
#' @examples
#' f <- c("bowtiePE.txt", "bowtieSE.txt")
#' bowtieLogs <- system.file("extdata", f, package = "ngsReports")
#' df <- importNgsLogs(bowtieLogs, type = "bowtie")
#'
#' @importFrom lubridate dminutes dhours dseconds parse_date_time
#'
#' @export
importNgsLogs <- function(x, type, which = 1) {

    x <- unique(x) # Remove any  duplicates
    stopifnot(file.exists(x)) # Check they all exist

    ## Check for a valid filetype
    possTypes <- c("bowtie", "bowtie2", "hisat2", "star", "duplicationMetrics")
    type <- match.arg(type, possTypes)
    ## Change to title case for easier parsing below
    type <- stringr::str_split_fixed(type, pattern = "", n = nchar(type))
    type[1] <- stringr::str_to_upper(type[1])
    type <- paste(type, collapse = "")

    ## Load the data
    data <- suppressWarnings(lapply(x, readLines))
    names(data) <- basename(x)

    ## Check validity
    vFun <- paste0(".isValid", type, "Log")
    validLogs <- vapply(data, eval(parse(text = vFun)), logical(1))
    if (any(!validLogs)) {
        failed <- names(validLogs)[!validLogs]
        stop(paste("Incorrect file structure for:", failed , collapse = "\n"))
    }

    ## Parse the data
    if (is.character(which)) which <- paste0("'", which, "'")
    pFun <- paste0(".parse", type, "Logs(data, which = ", which, ")")
    df <- eval(parse(text = pFun))

    ## Return the tibble
    as_tibble(df)

}

#' @title Check for correct structure of supplied Bowtie log files
#' @description Check for correct structure of supplied Bowtie log files after
#' reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.isValidBowtieLog <- function(x){
    nLines <- length(x)
    fields <- c(
        "Time loading forward index",
        "Time loading mirror index:",
        "Seeded quality full-index search:",
        "# reads processed:",
        "# reads with at least one reported alignment:",
        "# reads that failed to align:",
        "Time searching:",
        "Overall time:"
    )
    chk <- vapply(fields, function(f){any(grepl(f, x))}, logical(1))
    all(chk)
}

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
.isValidBowtie2Log <- .isValidHisat2Log

#' @title Check for a valid Star Alignment log
#' @description Checks internal structure of the parsed log file
#' @param x Character vector containing parsed log file using the function
#' readLines
#' @return logical(1)
#' @keywords internal
.isValidStarLog <- function(x){
    ## Just check for the key fields
    chkStart <- grepl("Started job on", x[1])
    chkUniq <- any(grepl("UNIQUE READS:", x))
    chkMulti <- any(grepl("MULTI-MAPPING READS:", x))
    chkUnmapped <- any(grepl("UNMAPPED READS:", x))
    all(chkStart, chkUniq, chkMulti, chkUnmapped)
}

#' @title Check for a valid Duplication Metrics log
#' @description Checks internal structure of the parsed log file
#' @param x Character vector containing parsed log file using the function
#' readLines
#' @return logical(1)
#' @keywords internal
.isValidDuplicationMetricsLog <- function(x){

    ## Check the METRICS CLASS data
    metricsHeader <- grep("METRICS CLASS\tpicard.sam.DuplicationMetrics", x)
    stopifnot(length(metricsHeader) == 1) # Must appear only once
    metricCols <- c("LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED",
                    "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS",
                    "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
                    "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION",
                    "ESTIMATED_LIBRARY_SIZE")
    ## Check the column names in the log match the expected columns
    checkMetCols <- all(names(.splitByTab(x[metricsHeader + 1])) == metricCols)

    ## Check the HISTOGRAM data
    histHeader <- grep("HISTOGRAM\tjava.lang.Double", x)
    stopifnot(length(histHeader) == 1) # Must appear only once
    histCols <- c("BIN", "VALUE")
    checkHistCols <- all(names(.splitByTab(x[histHeader + 1])) == histCols)

    all(checkMetCols, checkHistCols)
}

#' @title Parse data from Bowtie log files
#' @description Parse data from Bowtie log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseBowtieLogs <- function(data, ...){

    ## This will have already been through the validity checks
    df <- lapply(data, function(x){
        x <- gsub("# ", "", x)
        ## Treat the line containing the total reads separately
        total <- grep("Reported .* alignments", x = x, value = TRUE)
        x <- setdiff(x, total)
        ## Split the remaining lines into a nx2 matrix,
        ## then change to titles and remove spaces/dashes
        x <- stringr::str_split_fixed(x, pattern = ": ", 2)
        x[,1] <- stringr::str_to_title(x[,1])
        x[,1] <- gsub("( |-)", "_", x[,1])
        ## Remove those percentages from some of the fields
        x[,2] <- gsub("(.+) \\(.+\\)", "\\1", x[,2])
        ## Return a data.frame
        df <- structure(as.list(x[,2]), names = x[,1])
        df <- as.data.frame(df, stringsAsFactors = FALSE)
    })
    df <- dplyr::bind_rows(df)

    ## Some colnames may have a flag from the original bowtie code
    ## This line will correct that after the title case conversion above
    names(df) <- gsub("__(.)", "_-\\L\\1", names(df), perl = TRUE)

    ## Reformat the columns as integers and durations
    timeCols <- grepl("(Time|Full_Index_Search)", names(df))
    df[!timeCols] <- suppressWarnings(lapply(df[!timeCols], as.integer))
    df[timeCols] <- lapply(df[timeCols], function(x){
        x <- stringr::str_split_fixed(x, ":", 3)
        x <- as.numeric(x)
        x <- matrix(x, ncol = 3)
        dhours(x[,1]) + dminutes(x[,2]) + dseconds(x[,3])
    })
    ## Add the filename, reorder the columns & return a tibble
    df$Filename <- names(data)
    df <- dplyr::select(
        df,
        "Filename",
        tidyselect::contains("Reads"),
        tidyselect::contains("Time"),
        tidyselect::everything()
    )

    df
}

#' @title Parse data from HISAT2 log files
#' @description Parse data from HISAT2 log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseHisat2Logs <- function(data, ...){

    df <- lapply(data, function(x){

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

        out <- list(
            Total_Reads = totReads,
            Not_Aligned = noAln,
            Unique_Unpaired = uniq,
            Multiple_Unpaired = mult,
            Paired_Reads = pairReads,
            Unique_In_Pairs = uniqPairs,
            Multiple_In_Pairs = multPairs,
            Unique_Discordant_Pairs = uniqDiscord
        )
        ## Set all values as integers
        out <- lapply(out, as.integer)
        as.data.frame(out)
    })

    ## Bind all reults together
    df <- bind_rows(df)

    ## Add a final column
    df$Alignment_Rate <-
        1 - df$Not_Aligned / (df$Total_Reads + df$Paired_Reads)
    df$Filename <- names(data)
    dplyr::select(
        df,
        "Filename",
        tidyselect::ends_with("Reads"),
        tidyselect::contains("Unique"),
        tidyselect::contains("Multiple"),
        tidyselect::everything()
    )

}
.parseBowtie2Logs <- .parseHisat2Logs

#' @title Parse data from STAR log files
#' @description Parse data from STAR log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return tibble
#' @keywords internal
.parseStarLogs <- function(data, ...){
    ## Reformat as a data.frame / tibble
    df <- lapply(data, function(x){
        ## Split on '|\t'
        x <- stringr::str_split_fixed(x, pattern = "\\|\t", n = 2)
        ## Remove whitespace
        x <- apply(x, MARGIN = 2, stringr::str_trim)
        ## Remove blank rows
        x <- x[rowSums(x == "") == 0,]
        ## Tidy up the percent signs & brackets
        x[, 1] <- gsub("%", "Percent", x[,1])
        x <- apply(x, 2, stringr::str_remove_all, "[%\\(\\)]")
        ## Define the column names
        cols <- stringr::str_to_title(x[,1])
        cols <- gsub("[ :,]+", "_", cols)
        cols <-
            gsub("([ACGTacgt]{2}/[ACGTacgt]{2})", "\\U\\1", cols, perl = TRUE)
        ## Form a list then tibble
        out <- as.list(x[,2])
        names(out) <- cols
        ## Using tibble here preserves column names in a nicer format
        as_tibble(out)
    })
    df <- dplyr::bind_rows(df)

    ## Reformat the time columns
    timeCols <- grepl("On$", names(df))
    df[timeCols] <-
        lapply(df[timeCols], parse_date_time, orders = "b! d! HMS")

    ## Reformat the integer columns
    intCols <- grepl("Number", names(df))
    df[intCols] <- lapply(df[intCols], as.integer)

    ## Set the remaining columns as numeric
    df[!intCols & !timeCols] <- lapply(df[!intCols & !timeCols], as.numeric)

    ## Add the filename & additional columns
    df$Filename <- names(data)
    df$Mapping_Duration <- df$Finished_On - df$Started_Mapping_On
    totMapped <-
        df$Uniquely_Mapped_Reads_Number +
        df$Number_Of_Reads_Mapped_To_Multiple_Loci
    df$Total_Mapped_Percent <- 100*totMapped / df$Number_Of_Input_Reads

    ## Return the output
    dplyr::select(
        df,
        "Filename",
        tidyselect::starts_with("Total"),
        tidyselect::contains("Input"),
        tidyselect::contains("Mapped"),
        tidyselect::contains("Splice"),
        tidyselect::ends_with("On"),
        tidyselect::contains("Mapping"),
        tidyselect::everything()
    )

}

#' @title Parse data from STAR log files
#' @description Parse data from STAR log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param which which element of the log file to return. Can be 1:2, "metrics"
#' or "histogram"
#' @return tibble
#' @keywords internal
.parseDuplicationMetricsLogs <- function(data, which = 1){

    stopifnot(which %in% c(1, 2, "metrics", "histogram"))

    ## Collect the metrics from all files as a tibble
    metrics <- lapply(data, function(x){
        ## Find the library name
        libName <- grep(
            "picard.sam.markduplicates.MarkDuplicates.+INPUT=", x, value = TRUE
        )
        libName <- gsub(".+INPUT=\\[(.+)\\] OUTPUT.+", "\\1", libName[[1]])
        ## Find the header. The next two rows will be the colnames + data
        metHeader <- grep("METRICS CLASS\tpicard.sam.DuplicationMetrics", x)
        df <- .splitByTab(x[seq(metHeader + 1, by = 1, length.out = 2)])
        df$LIBRARY <- basename(libName)
        df
    })
    metrics <- dplyr::bind_rows(metrics)
    ## Now ensure the correct types
    metrics$PERCENT_DUPLICATION <- as.numeric(metrics$PERCENT_DUPLICATION)
    intCols <- setdiff(colnames(metrics), c("LIBRARY", "PERCENT_DUPLICATION"))
    metrics[intCols] <- lapply(metrics[intCols], as.integer)
    metrics <- as_tibble(metrics)

    ## Collect the histogram data from all files as a tibble
    histData <- lapply(data, function(x){
        ## Find the library name
        libName <- grep(
            "picard.sam.markduplicates.MarkDuplicates.+INPUT=", x, value = TRUE
        )
        libName <- gsub(".+INPUT=\\[(.+)\\] OUTPUT.+", "\\1", libName[[1]])

        ## Find the header then remove up until that line
        histHeader <- grep("HISTOGRAM\tjava.lang.Double", x)
        x <- x[-seq_len(histHeader)]
        ## Remove any blank lines (this is the last line)
        x <- x[!grepl("^$", x)]
        df <- .splitByTab(x)
        df$LIBRARY <- basename(libName)
        dplyr::select(df, "LIBRARY", tidyselect::everything())
    })
    histData <- dplyr::bind_rows(histData)
    ## Ensure the correct types
    histData$BIN <- as.integer(histData$BIN)
    histData$VALUE <- as.numeric(histData$VALUE)
    histData <- as_tibble(histData)

    ## Setup the output, then format the column names
    out <- list(metrics = metrics, histogram = histData)
    lapply(out, function(x){
        colnames(x) <- stringr::str_replace_all(colnames(x), "_", " ")
        colnames(x) <- stringr::str_to_title(colnames(x))
        x
    })
    out[[which]]
}
