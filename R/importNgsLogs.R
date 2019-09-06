#' @title Import Various NGS-related log files
#'
#' @description Imports NGS-related log files such as those generated from
#' stderr
#'
#' @details Imports one or more log files as output by tools such as:
#' \code{bowtie}, \code{bowtie2}, \code{featureCounts}, \code{Hisat2},
#' \code{STAR}, \code{picard MarkDuplicates}, \code{cutadapt},
#' \code{Adapter Removal}, \code{trimmomatic} \code{quast} or \code{busco}.
#' \code{autoDetect} can be used to detect the log type by parsing the file.
#'
#' The featureCounts log file corresponds to the \code{counts.out.summary},
#' not the main \code{counts.out} file.
#'
#' Whilst most log files return a single tibble, some are more complex
#' with multiple modules.
#'
#' \code{adapterRemoval} can return one of four modules (which = 1:4),.
#' When calling by name, the possible values are sequences, settings,
#' statistics or distribution.
#' Partial matching is implemented.
#'
#' \code{cutadapt} can return one of five modules (which = 1:5).
#' When calling by name the possible modules are summary, adapter1, adapter2,
#' adapter3 or overview.
#' Note that adapter2/3 may be missing from these files depending on the nature
#' of your data.
#' If cutadapt log files are obtained using \code{report=minimal}, all supplied
#' log files must be of this format and no modules can be returned.
#'
#' \code{duplicationMetrics} will return either the metrics of histogram.
#' These can be requested by setting which as 1 or 2, or naming either module.
#'
#' @param x \code{character}. Vector of filenames. All log files must be of the
#' same type. Duplicate file paths will be silently ignored.
#' @param type \code{character}. The type of file being imported. Can be one of
#' \code{bowtie}, \code{bowtie2}, \code{hisat2}, \code{star},
#' \code{featureCounts}, \code{duplicationMetrics}, \code{cutadapt},
#' \code{adapterRemoval}, \code{quast} or \code{busco}.
#' Defaults to \code{type = "auto"} which will automatically detect the file
#' type for all implemented types.
#' @param which Which element of the parsed object to return. Ignored in all
#' file types except when \code{type} is set to duplicationMetrics, cutadapt or
#' adapterRemoval. See details for possible values
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
#' @importFrom tidyselect contains everything starts_with ends_with
#'
#' @export
importNgsLogs <- function(x, type = "auto", which) {

    x <- unique(x) # Remove any  duplicates
    stopifnot(length(x) > 0) # Fail on empty vector
    stopifnot(file.exists(x)) # Check they all exist

    ## Check for a valid filetype
    possTypes <- c(
        "adapterRemoval",
        "bowtie",
        "bowtie2",
        "busco",
        "cutadapt",
        "duplicationMetrics",
        "featureCounts",
        "hisat2",
        "quast",
        "samtoolsStats",
        "star",
        "trimmomatic"
    )

    type <- match.arg(type, c("autoDetect", possTypes))


    ## Sort out the which argument
    if (missing(which)) {
        which <- 1L
        if (type == "adapterRemoval") {
            message("Defaulting to the statistics module for Adapter Removal")
            which <- 3L
        }
    }

    ## Load the data
    data <- suppressWarnings(lapply(x, readLines))
    names(data) <- basename(x)

    ## adding auto detect
    if (type == "autoDetect") {
        type <- unlist(lapply(data, .getToolName, possTypes = possTypes))
        type <- unique(type)
    }

    ## Change to title case for easier parsing below
    type <- stringr::str_split_fixed(type, pattern = "", n = nchar(type))
    type[1] <- stringr::str_to_upper(type[1])
    type <- paste(type, collapse = "")

    ## Check validity using the helpers defined as private functions
    vFun <- paste0(".isValid", type, "Log")
    validLogs <- vapply(data, eval(parse(text = vFun)), logical(1))
    if (any(!validLogs)) {
        failed <- names(validLogs)[!validLogs]
        stop(paste("Incorrect file structure for:", failed , collapse = "\n"))
    }

    ## The two types of cutadapt output also need to be checked for consistency
    if (type == "Cutadapt") {
        ln <- vapply(data, length, integer(1))
        nMinimal <- sum(ln == 2) # nLines for a minimal report
        if (!nMinimal %in% c(0, length(data))) stop(
            paste0(
                "Some, but not all, reports are minimal format.\n",
                "Please load different format reports separately"
            )
        )

    }

    ## Parse the data using the helpers defined as private functions
    if (is.character(which)) which <- paste0("'", which, "'")
    pFun <- paste0(".parse", type, "Logs(data, which = ", which, ")")
    df <- eval(parse(text = pFun))

    ## Return a tibble
    as_tibble(df)

}



#' @title Identify tool name
#' @description Identify tool name for log files after
#' reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.getToolName <- function(x, possTypes){

    logiList <- lapply(possTypes, function(types){
        ## Change to title case for easier parsing below
        types <- stringr::str_split_fixed(types, pattern = "", n = nchar(types))
        types[1] <- stringr::str_to_upper(types[1])
        types <- paste(types, collapse = "")
        vFun <- paste0(".isValid", types, "Log(x)")
        validLogs <- eval(parse(text = paste0(vFun)))

    })
    logi <- unlist(logiList)

    type <- possTypes[logi]
    ## added for autodetect purposes
    if (all(c("hisat2", "bowtie2") %in% type)) type <- "bowtie2"
    type
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

    if (length(metricsHeader) == 1) { # Must appear only once
        metricCols <- c("LIBRARY", "UNPAIRED_READS_EXAMINED",
                        "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS",
                        "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES",
                        "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES",
                        "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE")
        ## Check the column names in the log match the expected columns
        checkMetCols <-
            all(names(.splitByTab(x[metricsHeader + 1])) == metricCols)

        ## Check the HISTOGRAM data
        histHeader <- grep("HISTOGRAM\tjava.lang.Double", x)
        stopifnot(length(histHeader) == 1) # Must appear only once
        histCols <- c("BIN", "VALUE")
        checkHistCols <- all(names(.splitByTab(x[histHeader + 1])) == histCols)

        all(checkMetCols, checkHistCols)

    } else FALSE #give false value if missing

}

#' @title Check for a valid AdapterRemoval log
#' @description Checks internal structure of the parsed log file
#' @param x Character vector containing parsed log file using the function
#' readLines
#' @return logical(1)
#' @keywords internal
.isValidAdapterRemovalLog <- function(x){

    ## Check the first line begins with AdapterRemoval
    checkAR <- grepl("AdapterRemoval", x[[1]])

    ## This should contain four modules, all of which are contained
    ## between square brackets
    modNames <- c(
        "[Adapter sequences]",
        "[Adapter trimming]",
        "[Trimming statistics]",
        "[Length distribution]"
    )
    checkModNames <- stringr::str_subset(x[[1]], "^\\[.+") == modNames

    all(checkAR, checkModNames)
}

#' @title Check for a valid cutadapt log
#' @description Checks internal structure of the parsed log file
#' @param x Character vector containing parsed log file using the function
#' readLines
#' @return logical(1)
#' @keywords internal
.isValidCutadaptLog <- function(x){

    ## These can take two forms, minimal & full
    ## The minimal report is like a tab-delimited data.frame
    ## whilst the full report has a more complex structure
    isMinimal <- length(x) == 2
    reqCols <- c(
        "status", "in_reads", "in_bp", "too_short", "too_long",
        "too_many_n",  "out_reads", "w/adapters", "qualtrim_bp", "out_bp"
    )
    if (isMinimal) {
        ## As this is essentially a df, just check the columns
        df <- tryCatch(.splitByTab(x))
        chkCols <- all(reqCols %in% colnames(df))
        return(TRUE)
    }
    ## Check the complete file
    chkCutAdapt <- grepl("cutadapt", x[[1]])
    chkCols <- list(
        status = any(grepl("Finished in", x)),
        in_reads = any(grepl("Total reads processed", x)),
        in_bp = any(grepl("Total basepairs processed", x)),
        too_short = any(grepl("Reads that were too short", x)),
        ## Too long may not be in the file
        too_many_n = any(grepl("Reads with too many N", x)),
        out_reads = any(grepl("Reads written \\(passing filters\\)", x)),
        `w/adapters` = any(grepl("Reads with adapters", x)),
        qualtrim_bp = any(grepl("Quality-trimmed", x)),
        out_bp = any(grepl("Total written \\(filtered\\)", x))
    )
    chkCols <- unlist(chkCols)
    all(c(chkCutAdapt, chkCols))

}

#' @title Check for a valid featureCounts Summary
#' @description Checks internal structure of the parsed file
#' @param x Character vector containing parsed log file using the function
#' readLines
#' @return logical(1)
#' @keywords internal
.isValidFeatureCountsLog <- function(x){

    if (all(grepl("\t", x))){ ### check all lines have a tab,
        ## required by .splitbyTab

        ## Check the Status column
        x <- .splitByTab(x)
        chkCol1 <- colnames(x)[[1]] == "Status"
        chkTypes <- FALSE
        if (chkCol1) {
            ## We can only do this if the Status column checks out
            vals <- c(
                "Assigned", "Unassigned_Ambiguity", "Unassigned_MultiMapping",
                "Unassigned_NoFeatures", "Unassigned_Unmapped",
                "Unassigned_MappingQuality",
                "Unassigned_FragmentLength", "Unassigned_Chimera",
                "Unassigned_Secondary", "Unassigned_Nonjunction",
                "Unassigned_Duplicate"
            )
            chkTypes <- all(vals %in% x[["Status"]])
        }

        all(chkTypes, chkCol1)
    }
    else FALSE ## print false if not

}

#' @title Check for correct structure of supplied BUSCO log files
#' @description Check for correct structure of supplied BUSCO log files after
#' reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.isValidBuscoLog <- function(x){

    fields <- c(
        "# BUSCO version is:",
        "# The lineage dataset is:",
        "# To reproduce this run:",
        "# Summarized benchmarking in BUSCO notation for file ",
        "# BUSCO was run in mode:"
    )
    chk <- vapply(fields, function(f){any(grepl(f, x))}, logical(1))
    all(chk)
}

#' @title Check for correct structure of supplied Trimmomatic log files
#' @description Check for correct structure of supplied Trimmomatic log files
#' after reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.isValidTrimmomaticLog <- function(x){
    n <- length(x)
    checkL1 <- grepl("Trimmomatic[PS]E: Started with arguments", x[[1]])
    checkMain <- grepl("Input.+ Surviving.+ Dropped.+", x[[n - 1]])
    checkLast <- grepl("Trimmomatic[PS]E: Completed successfully", x[[n]])
    all(checkL1, checkMain, checkLast)
}


#' @title Check for correct structure of supplied Quast log files
#' @description Check for correct structure of supplied Quast log files after
#' reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.isValidQuastLog <- function(x){

    fields <- c(
        "# contigs",
        "Largest contig",
        "N50",
        "N75",
        "Total length",
        "L50",
        "L75",
        "# N's per 100 kbp"
    )
    chk <- vapply(fields, function(f){any(grepl(f, x))}, logical(1))
    all(chk)
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
    dplyr::select(
        df, "Filename", contains("Reads"), contains("Time"), everything()
    )

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
        ends_with("Reads"),
        contains("Unique"),
        contains("Multiple"),
        everything()
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
        starts_with("Total"),
        contains("Input"),
        contains("Mapped"),
        contains("Splice"),
        ends_with("On"),
        contains("Mapping"),
        everything()
    )

}

#' @title Parse data from Picard duplicationMetrics log files
#' @description Parse data from Picard duplicationMetrics log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param which which element of the log file to return.
#' Can be 1:2, "metrics" or "histogram"
#' @return tibble
#' @keywords internal
.parseDuplicationMetricsLogs <- function(data, which = 1){

    if (is.character(which))
        which <- match.arg(which, c("metrics", "histogram"))
    if (is.numeric(which)) stopifnot(which %in% c(1, 2))

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
        dplyr::select(df, "LIBRARY", everything())
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

#' @title Parse data from Adapter Removal log files
#' @description Parse data from Adapter Removal log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param which which element of the log file to return.
#' Can be 1:4, "sequences", "settings", "statistics" or "distribution"
#' @return tibble
#' @keywords internal
#' @importFrom stringr str_split_fixed
.parseAdapterRemovalLogs <- function(data, which = 3){

    ## These are renamed for easier understanding.
    ## The full mapping of titles are:
    ## 1: sequences = Adapter sequences
    ## 2: settings = Adapter trimming
    ## 3: statistics =  Trimming statistics
    ## 4: distribution = Length distribution
    ## No parsing of the version number has been implemented (1st row)
    ## Note that for each file, 'distribution' will return a tibble
    ## with nrow > 1
    modNames <- c("sequences", "settings", "statistics", "distribution")
    if (is.character(which)) which <- match.arg(which, modNames)
    if (is.numeric(which)) stopifnot(which %in% seq_len(4))

    .parseArSingle <- function(x, which){
        ## Remove the blank elements
        x <- setdiff(x, "")
        ## Find where the modules start
        vals2Mods <- cumsum(stringr::str_detect(x, "^\\["))

        ## Just do the modules manually
        ## Start with the sequences
        sequences <- str_split_fixed(
            string = x[vals2Mods == 1][-1],
            pattern = ": ",
            n =  2
        )
        sequences <- as_tibble(
            split(sequences[,2], f = sequences[,1])
        )
        ## Now the settings
        settings <- str_split_fixed(
            string = x[vals2Mods == 2][-1],
            pattern = ": ",
            n =  2
        )
        settings <- split(
            x = settings[,2],
            f = factor(settings[,1], levels = settings[,1])
        )
        ## Convert to the correct value type
        intCols <- grepl("(score |value|length|Minimum)", names(settings))
        intCols <- intCols & !grepl("Maximum", names(settings))
        settings[intCols] <- lapply(settings[intCols], as.integer)
        numCols <- grepl("threshold", names(settings))
        settings[numCols] <- lapply(settings[numCols], as.numeric)
        settings <- as_tibble(settings)
        settings[["RNG seed"]] <- ifelse(
            settings[["RNG seed"]] == "NA",
            NA_integer_,
            as.integer(settings[["RNG seed"]])
        )
        ## Statistics
        statistics <- str_split_fixed(x[vals2Mods == 3][-1], ": ", 2)
        statistics <- split(
            x = statistics[,2],
            f = factor(statistics[,1], levels = statistics[,1])
        )
        statistics <- lapply(statistics, as.numeric)
        statistics <- as_tibble(statistics)
        ## Distribution
        distribution <- .splitByTab(x[vals2Mods == 4][-1])
        distribution <- lapply(distribution, as.integer)
        distribution <- as_tibble(distribution)

        ## Setup the output, then format the column names
        out <- list(
            sequences = sequences,
            settings = settings,
            statistics = statistics,
            distribution = distribution
        )
        out[[which]]
    }

    arOut <- lapply(data, .parseArSingle, which = which)
    arOut <- lapply(names(arOut), function(x){
        x <- dplyr::mutate(arOut[[x]], Filename = x)
        dplyr::select(x, Filename, everything())
    })
    ## Now bind all tibbles & returns
    bind_rows(arOut)

}

#' @title Parse data from cutadapt log files
#' @description Parse data from cutadapt log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param which which element of the log file to return.
#' Can be summary, adapter1, adapter2, adapter3 or overview, or any integrer in
#' 1:5
#' @return tibble
#' @keywords internal
#' @importFrom tidyselect everything
#' @importFrom dplyr mutate_at mutate_if
#' @importFrom stringr str_replace_all str_remove_all str_extract
.parseCutadaptLogs <- function(data, which = 1){

    ## Possible modules for the full version are:
    ## 1: summary (this is identical to minimal)
    ## 2: adapter1/2/3
    ## 5: overview
    ## No parsing of the version number has been implemented (1st row)
    modNames <- c(
        "summary", "adapter1", "adapter2", "adapter3", "overview"
    )
    if (is.character(which)) which <- match.arg(which, modNames)
    if (is.numeric(which)) {
        stopifnot(which %in% seq_len(5))
        which <- modNames[which]
    }

    .parseCutAdaptSingle <- function(x, which){

        ## Check if it is a minimal report & return the df if TRUE
        isMinimal <- length(x) == 2
        if (isMinimal) {
            df <- .splitByTab(x)
            ## Values may sometimes be > 2^31 so convert to doubles
            df <- dplyr::mutate_if(df, colnames(df) != "status", as.numeric)
            if (which != 1) message(
                paste0(
                    "Minimal report provided.",
                    "The 'which' argument will be ignored"
                )
            )
            return(df)
        }

        ## Otherwise parse the modules:
        ## Remove the blank elements
        x <- setdiff(x, "")
        ## Find where the modules start
        mods <- grepl("(=== Summary|=== Adapter|Overview)", x)
        foundMods <- tolower(str_remove_all(x[mods],"[= ]"))
        foundMods <- gsub("overview.+", "overview", foundMods)
        ## Check the requested module exists as there is some variability
        modExists <- which %in% foundMods
        if (!modExists) stop("The requested module (", which, ") is missing")

        ## Now split into the modules
        x <- split(x, f = cumsum(mods))
        ## Remove the first value from each
        x <- lapply(x, function(x){x[-1]})
        names(x) <- c("header", foundMods)

        ## Just do the modules manually
        out <- vector("list", length(x) - 1)
        names(out) <- names(x)[-1]

        ## Start with the summary (i.e. minimal)
        vals <- str_replace_all(x$summary, ".+ +([0-9,]+).+", "\\1")
        vals <- str_remove_all(vals, ",")
        vals <- as.numeric(vals)
        names(vals) <- str_replace_all(x$summary, "(.+):.+", "\\1")
        reqCols <- c(
            "in_reads", "in_bp", "too_short", "too_long",
            "too_many_n",  "out_reads", "w/adapters", "qualtrim_bp", "out_bp"
        )
        fullText <- c(
            "Total reads processed", "Total basepairs processed",
            "Reads that were too short", "Reads that were too long",
            "Reads with too many N", "Reads written (passing filters)",
            "Reads with adapters", "Quality-trimmed",
            "Total written (filtered)"
        )
        names(vals) <- reqCols[match(names(vals), fullText)]
        out[["summary"]] <-
            as_tibble(as.list(vals))[intersect(reqCols, names(vals))]

        ## Now the adapters
        out[grep("adapter", names(out))] <- lapply(
            x[grepl("adapter", names(x))],
            function(a){
                tibble(
                    sequence = str_extract(
                        a[grepl("Sequence", a)], "[ACGTN]+"
                    ),
                    type = str_replace_all(
                        a[grepl("Type", a)],
                        ".+Type: (.+); Length.+",
                        "\\1"
                    ),
                    length = as.integer(
                        str_replace_all(
                            a[grepl("Length", a)],
                            ".+Length: (.+); Trimmed.+",
                            "\\1"
                        )
                    ),
                    A = as.numeric(
                        str_extract(a[grepl("A:", a)], "[0-9\\.]+")
                    ) / 100,
                    C = as.numeric(
                        str_extract(a[grepl("C:", a)], "[0-9\\.]+")
                    ) / 100,
                    G = as.numeric(
                        str_extract(a[grepl("G:", a)], "[0-9\\.]+")
                    ) / 100,
                    `T` = as.numeric(
                        str_extract(a[grepl("T:", a)], "[0-9\\.]+")
                    ) / 100
                )
            }
        )

        ## The overview
        out[["overview"]] <- .splitByTab(x$overview)
        out[["overview"]] <- dplyr::mutate_if(
            out[["overview"]],
            vapply(
                out[["overview"]],
                function(x){suppressWarnings(all(!is.na(as.numeric(x))))},
                logical(1)
            ),
            as.numeric
        )
        out[["overview"]] <- as_tibble(out[["overview"]])

        ## Return the required module
        out[[which]]
    }

    caOut <- lapply(data, .parseCutAdaptSingle, which = which)
    caOut <- lapply(names(caOut), function(x){
        x <- dplyr::mutate(caOut[[x]], Filename = x)
        dplyr::select(x, Filename, everything())
    })
    ## Now bind all tibbles & returns
    bind_rows(caOut)

}

#' @title Parse data from featureCounts summary files
#' @description Parse data from featureCounts summary files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' Can be 1:4, "sequences", "settings", "statistics" or "distribution"
#' @return tibble
#' @keywords internal
.parseFeatureCountsLogs <- function(data, ...){

    out <- lapply(data, function(x){
        x <- .splitByTab(x)
        Status <- c() # Avoiding an R CMD check error
        x <- tidyr::gather(x, "Sample", "Total", -Status)
        x$Sample <- basename(x$Sample)
        x$Total <- as.integer(x$Total)
        x$Status <- factor(x$Status, levels = unique(x$Status))
        tidyr::spread(x, "Status", "Total")
    })

    out <- lapply(names(data), function(x){
        out[[x]]$Filename = x
        dplyr::select(out[[x]], Filename, everything())
    })
    bind_rows(out)
}


#' @title Parse data from BUSCO log files
#' @description Parse data from BUSCO log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseBuscoLogs <- function(data, ...){

    ## This will have already been through the validity checks
    df <- lapply(data, function(x){
        x <- gsub("# ", "", x)
        ## get numbers of each feild from lines
        single <-
            grep("Complete and single-copy BUSCOs \\(S\\)", x = x, value = TRUE)
        single <- as.integer(gsub(".*\t(.+)\tComp.*", "\\1", single))

        duplicated <-
            grep("Complete and duplicated BUSCOs \\(D\\)", x = x, value = TRUE)
        duplicated <- as.integer(gsub(".*\t(.+)\tComp.*", "\\1", duplicated))

        fragmented <- grep("Fragmented BUSCOs \\(F\\)", x = x, value = TRUE)
        fragmented <- as.integer(gsub(".*\t(.+)\tFrag.*", "\\1", fragmented))

        missing <- grep("Missing BUSCOs \\(M\\)", x = x, value = TRUE)
        missing <- as.integer(gsub(".*\t(.+)\tMiss.*", "\\1", missing))

        total <- sum(single, duplicated, fragmented, missing)

        name <- grep(
            "Summarized benchmarking in BUSCO notation for file ",
            x = x,
            value = TRUE
        )
        name <- gsub(".*file ", "\\1", name)
        df <- tibble(
            name = name,
            completeSingleCopy = single,
            completeDuplicated = duplicated,
            fragmented = fragmented,
            missing = missing
        )
    })
    df <- dplyr::bind_rows(df)

}

#' @title Parse data from trimmomatic log files
#' @description Parse data from trimmomatic log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... not used
#' @return tibble
#' @keywords internal
#' @importFrom tidyselect everything starts_with contains
.parseTrimmomaticLogs <- function(data, ...){

    .parseTrimmoSingle <- function(x){

        ## Initialise values which may or may not be present
        Input_Reads <- Input_Read_Pairs <- Surviving <- Both_Surviving <- NA
        Forward_Only_Surviving <- Reverse_Only_Surviving <- NA

        readType <- gsub(".+(PE|SE).+", "\\1", x[[1]])
        Illumina_Clip <- ifelse(
            grepl("ILLUMINACLIP", x[[2]]),
            gsub(".+ILLUMINACLIP:([^ ]+).+", "\\1", x[[2]]),
            NA
        )
        Leading <- ifelse(
            grepl("LEADING", x[[2]]),
            as.integer(gsub(".+LEADING:([0-9:]+).+", "\\1", x[[2]])),
            NA_integer_
        )
        Trailing <- ifelse(
            grepl("TRAILING", x[[2]]),
            as.integer(gsub(".+TRAILING:([0-9:]+).+", "\\1", x[[2]])),
            NA_integer_
        )
        Crop <- ifelse(
            grepl(" CROP", x[[2]]),
            as.integer(gsub(".+ CROP:([0-9:]+).+", "\\1", x[[2]])),
            NA_integer_
        )
        Head_Crop <- ifelse(
            grepl("HEADCROP", x[[2]]),
            as.integer(gsub(".+HEADCROP:([0-9:]+).+", "\\1", x[[2]])),
            NA_integer_
        )
        Sliding_Window <- ifelse(
            grepl("SLIDINGWINDOW", x[[2]]),
            gsub(".+SLIDINGWINDOW:([0-9:]+).+", "\\1", x[[2]]),
            NA
        )
        Min_Len <- ifelse(
            grepl("MINLEN", x[[2]]),
            gsub(".+MINLEN:([0-9:]+).+", "\\1", x[[2]]),
            NA
        )
        Max_Info <- ifelse(
            grepl("MAXINFO", x[[2]]),
            gsub(".+MAXINFO:([^ ]+).+", "\\1", x[[2]]),
            NA
        )
        Avg_Qual <- ifelse(
            grepl("AVGQUAL", x[[2]]),
            as.integer(gsub(".*AVGQUAL:([0-9]+).*", "\\1", x[[2]])),
            NA_integer_
        )
        Quality_Encoding <- grep("Quality encoding", x, value = TRUE)
        Quality_Encoding <-
            gsub("Quality encoding detected as ", "", Quality_Encoding)

        ## Get the line with the summary values
        valLine <- x[[length(x) - 1]]

        if (readType == "SE") {

            Input_Reads <- gsub("Input Reads: ([0-9]+).+", "\\1", valLine)
            Input_Reads <- as.integer(Input_Reads)

            Surviving <- gsub(".+Surviving: ([0-9]+).+", "\\1", valLine)
            Surviving <- as.integer(Surviving)

        }
        if (readType == "PE") {

            Input_Read_Pairs <-
                gsub("Input Read Pairs: ([0-9]+).+", "\\1", valLine)
            Input_Read_Pairs <- as.integer(Input_Read_Pairs)

            Both_Surviving <-
                gsub(".+Both Surviving: ([0-9]+).+", "\\1", valLine)
            Both_Surviving <- as.integer(Both_Surviving)

            Forward_Only_Surviving <-
                gsub(".+Forward Only Surviving: ([0-9]+).+", "\\1", valLine)
            Forward_Only_Surviving <- as.integer(Forward_Only_Surviving)

            Reverse_Only_Surviving <-
                gsub(".+Reverse Only Surviving: ([0-9]+).+", "\\1", valLine)
            Reverse_Only_Surviving <- as.integer(Reverse_Only_Surviving)
        }
        Dropped <- gsub(".+Dropped: ([0-9]+).+", "\\1", valLine)
        Dropped <- as.integer(Dropped)

        tibble(
            Type = readType,
            Input_Reads,
            Input_Read_Pairs,
            Surviving,
            Both_Surviving,
            Forward_Only_Surviving,
            Reverse_Only_Surviving,
            Dropped,
            Illumina_Clip,
            Sliding_Window,
            Max_Info,
            Leading,
            Trailing,
            Crop,
            Head_Crop,
            Min_Len,
            Avg_Qual,
            Quality_Encoding
        )

    }

    out <- lapply(data, .parseTrimmoSingle)
    out <- dplyr::bind_rows(out)
    out$Filename <- names(data)

    ## Many of the above values may be missing.
    ## Remove them if so using a quick tidy
    value <- c() # Avoiding an R CMD check NOTE
    out <- tidyr::gather(out, "key", "value", -1)
    out <- dplyr::filter(out, !is.na(value))
    out <- tidyr::spread(out, "key", "value")

    ## Return the final output
    dplyr::select(
        out,
        "Filename",
        "Type",
        starts_with("Input"),
        contains("Surviving"),
        everything()
    )

}


#' @title Parse data from BUSCO log files
#' @description Parse data from BUSCO log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseQuastLogs <- function(data, ...){

    ## This will have already been through the validity checks
    df <- lapply(data, function(x){
        x <- gsub("# ", "", x)
        ## get Assembly stats from lines
        totalL <- grep("Total length\t", x = x, value = TRUE)
        totalL <- as.integer(gsub(".*\t(.+)", "\\1", totalL))

        N50 <- grep("N50\t", x = x, value = TRUE)
        N50 <- as.integer(gsub(".*\t(.+)", "\\1", N50))

        N75 <- grep("N75\t", x = x, value = TRUE)
        N75 <- as.integer(gsub(".*\t(.+)", "\\1", N75))

        L50 <- grep("L50\t", x = x, value = TRUE)
        L50 <- as.integer(gsub(".*\t(.+)", "\\1", L50))

        L75 <- grep("L75\t", x = x, value = TRUE)
        L75 <- as.integer(gsub(".*\t(.+)", "\\1", L75))

        longest <- grep("Largest contig\t", x = x, value = TRUE)
        longest <- as.integer(gsub(".*\t(.+)", "\\1", longest))

        df <- tibble(
            totalLength = totalL,
            longestTig = longest,
            N50 = N50,
            N75 = N75,
            L50 = L50,
            L75 = L75
        )
    })
    df <- dplyr::bind_rows(df)

    dplyr::bind_cols(tibble(fileNames = names(data)), df)

}