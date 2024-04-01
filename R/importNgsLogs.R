#' @title Import Various NGS-related log files
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#' Imports NGS-related log files such as those generated from
#' stderr.
#'
#' @details Imports one or more log files as output by tools such as:
#' `bowtie`, `bowtie2`, `featureCounts`, `Hisat2`, `STAR`, `salmon`
#' `picard MarkDuplicates`, `cutadapt`, `flagstat`, `macs2Callpeak`,
#' `Adapter Removal`, `trimmomatic`, `rnaseqcMetrics`, `quast` or `busco`.
#' `autoDetect` can be used to detect the log type by parsing the file.
#'
#' The featureCounts log file corresponds to the `counts.out.summary`,
#' not the main `counts.out` file.
#'
#' Whilst most log files return a single tibble, some are more complex
#' with multiple modules.
#'
#' `adapterRemoval` can return one of four modules (which = 1:4),.
#' When calling by name, the possible values are sequences, settings,
#' statistics or distribution.
#' Partial matching is implemented.
#'
#' `cutadapt` can return one of five modules (which = 1:5).
#' When calling by name the possible modules are summary, adapter1, adapter2,
#' adapter3 or overview.
#' Note that adapter2/3 may be missing from these files depending on the nature
#' of your data.
#' If cutadapt log files are obtained using `report=minimal`, all supplied
#' log files must be of this format and no modules can be returned.
#'
#' `duplicationMetrics` will return either the metrics of histogram.
#' These can be requested by setting which as 1 or 2, or naming either module.
#'
#' @param x `character`. Vector of filenames. All log files must be of the
#' same type. Duplicate file paths will be silently ignored.
#' @param type `character`. The type of file being imported. Can be one of
#' `bowtie`, `bowtie2`, `hisat2`, `star`, `flagstat`,
#' `featureCounts`, `duplicationMetrics`, `cutadapt`, `umitoolsDedup`,
#' `macs2Callpeak`, `adapterRemoval`, `rnaseqcMetrics`, `quast`,
#' `salmonLibFormatCounts`, `salmonMetaInfo` or `busco`.
#' Defaults to `type = "auto"` which will automatically detect the file
#' type for all implemented types.
#' @param which Which element of the parsed object to return. Ignored in all
#' file types except when `type` is set to duplicationMetrics, cutadapt or
#' adapterRemoval. See details for possible values.
#' To return all elements, set this value to 'all'
#' @param stripPaths logical(1). Remove paths from the Filename column
#'
#' @return A `tibble`.
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
#' @importFrom stringr str_replace_all str_trim str_remove_all str_extract
#' @importFrom tidyr pivot_wider
#'
#' @export
importNgsLogs <- function(x, type = "auto", which, stripPaths = TRUE) {

    x <- unique(x) # Remove any  duplicates
    stopifnot(length(x) > 0) # Fail on empty vector
    stopifnot(file.exists(x)) # Check they all exist

    ## Check for a valid filetype
    possTypes <- c(
        "adapterRemoval",
        "bowtie", "bowtie2",
        "busco",
        "cutadapt",
        "duplicationMetrics",
        "featureCounts",
        "flagstat",
        "hisat2",
        "macs2Callpeak",
        "rnaseqcMetrics",
        "quast",
        "salmonLibFormatCounts",
        "salmonMetaInfo",
        "star",
        "trimmomatic",
        "umitoolsDedup"
    )
    type <- match.arg(type, c("autoDetect", possTypes))

    ## Load the data
    data <- suppressWarnings(lapply(x, readLines))
    if (stripPaths) x <- basename(x)
    names(data) <- x

    ## adding auto detect
    if (type == "autoDetect") {
        type <- unlist(lapply(data, .getToolName, possTypes = possTypes))
        type <- unique(type)
    }

    ## Sort out the which argument
    if (missing(which)) {
        which <- 1L
        if (type == "adapterRemoval") {
            message("Defaulting to the statistics module for Adapter Removal")
            which <- 3L
        }
    }

    ## Change to title case for easier parsing below
    type <- stringr::str_split_fixed(type, pattern = "", n = nchar(type))
    type[1] <- stringr::str_to_upper(type[1])
    type <- paste(type, collapse = "")

    ## Check validity using the helpers defined as private functions
    vFun <- paste0(".isValid", type, "Log")
    validLogs <- vapply(data, eval(parse(text = vFun)), logical(1))
    if (any(!validLogs)) {
        failed <- paste(names(validLogs)[!validLogs], collapse = "\n")
        stop("Incorrect file structure for:\n", failed)
    }

    ## The two types of cutadapt output also need to be checked for consistency
    if (type == "Cutadapt") {
        ln <- vapply(data, length, integer(1))
        nMinimal <- sum(ln == 2) # nLines for a minimal report
        if (!nMinimal %in% c(0, length(data))) stop(
            "Some, but not all, reports are minimal format.\n",
            "Please load different format reports separately"
        )

    }

    ## Parse the data using the helpers defined as private functions
    if (is.character(which)) which <- paste0("'", which, "'")
    pFun <- paste0(".parse", type, "Logs(data, which = ", which, ")")
    out <- eval(parse(text = pFun))
    ## Return a tibble if only one module
    if (!grepl("all", which)) out <- as_tibble(out)
    out

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
    if (length(type) == 0) stop("No matching file type was found")
    type
}


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

.isValidHisat2Log <- function(x){
    n <- length(x)
    chkLen <- length(x) > 0
    firstLine <- grepl("reads; of these:$", x[1])
    lastLine <- grepl("overall alignment rate$", x[n])
    noAln <- sum(grepl("aligned.*0 times$", x)) %in% seq_len(2)
    alnExact <- sum(grepl("aligned.*exactly 1 time$", x)) %in% seq_len(2)
    alnG1 <- sum(grepl("aligned.*>1 times$", x)) %in% seq_len(2)
    all(c(chkLen, firstLine, lastLine, noAln, alnExact, alnG1))
}
.isValidBowtie2Log <- .isValidHisat2Log

.isValidStarLog <- function(x){
    ## Just check for the key fields
    chkStart <- grepl("Started job on", x[1])
    chkUniq <- any(grepl("UNIQUE READS:", x))
    chkMulti <- any(grepl("MULTI-MAPPING READS:", x))
    chkUnmapped <- any(grepl("UNMAPPED READS:", x))
    all(chkStart, chkUniq, chkMulti, chkUnmapped)
}

.isValidDuplicationMetricsLog <- function(x){

    ## Check the METRICS CLASS data
    metricsHeader <- grep("METRICS CLASS\t.+DuplicationMetrics", x)

    if (length(metricsHeader) == 1) { # Must appear only once
        metricCols <- c(
            "LIBRARY", "UNPAIRED_READS_EXAMINED",
            "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS",
            "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES",
            "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES",
            "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE"
        )
        ## Check the column names in the log match the expected columns
        checkMetCols <- all(
            names(.splitByTab(x[metricsHeader + 1])) == metricCols
        )

        ## Check the HISTOGRAM data
        checkHistCols <- TRUE
        if (any(grepl("HISTOGRAM", x))) {
            histHeader <- grep("HISTOGRAM\tjava.lang.Double", x)
            stopifnot(length(histHeader) == 1) # Must appear only once
            histCols <- c("BIN")
            nm <- names(.splitByTab(x[histHeader + 1]))[seq_along(histCols)]
            checkHistCols <- all(nm == histCols)
        }
        all(checkMetCols, checkHistCols)

    } else FALSE #give false value if missing

}

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
        in_reads = any(grepl("Total (reads|read pairs) processed", x)),
        in_bp = any(grepl("Total basepairs processed", x)),
        too_short = any(grepl("(Reads|Pairs) that were too short", x)),
        ## Too long may not be in the file
        too_many_n = any(grepl("(Reads|Pairs) with too many N", x)),
        out_reads = any(
            grepl("(Reads|Pairs) written \\(passing filters\\)", x)
        ),
        `w/adapters` = any(grepl("Read.+ with adapter", x)),
        out_bp = any(grepl("Total written \\(filtered\\)", x))
    )
    chkCols <- unlist(chkCols)
    all(c(chkCutAdapt, chkCols))

}

.isValidFeatureCountsLog <- function(x){

    if (all(grepl("\t", x))) { ### check all lines have a tab,
        ## required by .splitbyTab

        ## Check the Status column
        x <- .splitByTab(x)
        chkCol1 <- colnames(x)[[1]] == "Status"
        chkTypes <- FALSE
        if (chkCol1) {
            ## We can only do this if the Status column checks out
            vals <- c(
                "Assigned",
                "Unassigned_Ambiguity",
                "Unassigned_MultiMapping",
                "Unassigned_NoFeatures",
                "Unassigned_Unmapped",
                "Unassigned_MappingQuality",
                "Unassigned_FragmentLength",
                "Unassigned_Chimera",
                "Unassigned_Secondary",
                "Unassigned_Duplicate"
            )
            chkTypes <- all(vals %in% x[["Status"]])
        }

        all(chkTypes, chkCol1)
    }
    else FALSE ## print false if not

}

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

.isValidTrimmomaticLog <- function(x){
    n <- length(x)
    checkL1 <- grepl("Trimmomatic[PS]E: Started with arguments", x[[1]])
    checkMain <- grepl("Input.+ Surviving.+ Dropped.+", x[[n - 1]])
    checkLast <- grepl("Trimmomatic[PS]E: Completed successfully", x[[n]])
    all(checkL1, checkMain, checkLast)
}

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

.isValidFlagstatLog <- function(x){

    ## Every line must have a '+' symbol
    nLines <- length(x)
    allPlus <- sum(grepl("\\+", x)) == nLines
    ## Check for required text in the first line
    firstLine <- grepl("QC-passed reads", x[[1]])
    all(allPlus, firstLine)
}

.isValidMacs2CallpeakLog <- function(x){
    hasCmd <- grepl("callpeak", x[[2]])
    hasArgs <- any(grepl("ARGUMENTS LIST", x))
    hasBlankSep <- any(str_trim(x) == "")
    hasTagSize <- any(grepl("tag size", x))
    hasFragLength <- any(grepl("fragment length", x))
    all(hasCmd, hasArgs, hasBlankSep, hasTagSize, hasFragLength)
}

.isValidUmitoolsDedupLog <- function(x){
    isUmiTools <- grepl("UMI-tools", x[[1]])
    isDedup <- grepl("output generated by dedup", x[[2]])
    n <- length(x)
    isComplete <- grepl("job finished", x[[n]])
    hasData <- sum(grepl("INFO", x)) > 0
    all(isUmiTools, isDedup, isComplete, hasData)
}

.isValidRnaseqcMetricsLog <- function(x){
    ## Just check some key fields which should distinguish these files
    hasSampleFirst <- which(grepl("Sample", x)) == 1
    rt_pattern <- "(Exon|Intron|Intergen|Intragen)ic Rate"
    hasAllRates <- sum(grepl(paste0("^", rt_pattern), x)) == 4
    hasAllHQRates <- sum(grepl(paste("High Quality", rt_pattern), x)) == 4
    tot_pattern <- "Total (Bases|Reads|Mapped)"
    hasAllTotals <- sum(grepl(tot_pattern, x)) == 3
    all(hasSampleFirst, hasAllRates, hasAllHQRates, hasAllTotals)
}

.isValidSalmonLibFormatCountsLog <- function(x){
    ## Check a few key columns
    has_read_files <- any(grepl("read_files", x))
    has_format <- any(grepl("expected_format", x))
    has_comp_frags <- sum(grepl("compatible_frag", x)) == 2
    has_sfr <- sum(grepl("[MOI]S[FR]", x)) >= 6
    has_assigned <- any(grepl("assigned_frag", x))
    all(has_format, has_read_files, has_sfr, has_comp_frags, has_assigned)
}

.isValidSalmonMetaInfoLog <- function(x){
    ## Check a few key columns
    has_vers <- any(grepl("salmon_version", x))
    frag_len <- sum(grepl("frag_length", x)) == 2
    has_libtype <- any(grepl("library_types", x))
    has_bias <- sum(grepl("bias", x)) == 3
    has_hash <- sum(grepl("hash", x)) == 6
    has_time <- sum(grepl("time", x)) == 2
    all(has_vers, frag_len, has_libtype, has_bias, has_hash, has_time)
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
        noAln <- c(
            gsub("([0-9]*) .+ aligned 0 times", "\\1", noAln), NA_character_
        )[1]
        uniq <- x[grep("aligned exactly 1 time$", x)]
        uniq <- c(
            gsub("([0-9]*) .+ aligned exactly 1 time", "\\1", uniq),
            NA_character_
        )[1]
        mult <- x[grep("aligned >1 times", x)]
        mult <- c(
            gsub("([0-9]*) .+ aligned >1 times", "\\1", mult), NA_character_
        )[1]
        alnRate <- x[grep("overall alignment rate", x)]
        alnRate <- gsub("([0-9\\.]+)% overall.+", "\\1", alnRate)
        alnRate <- as.numeric(alnRate) / 100

        ## Get the paired only fields
        pairReads <- uniqPairs <- multPairs <- uniqDiscord <- NA_character_
        if (paired) {
            pairReads <- x[grep("were paired; of these:$", x)]
            pairReads <- c(
                gsub("([0-9]*) .+ of these:", "\\1", pairReads), NA_character_
            )[1]
            uniqPairs <- x[grep("aligned concordantly exactly 1 time$", x)]
            uniqPairs <- c(
                gsub("([0-9]*) .+ exactly 1 time", "\\1", uniqPairs),
                NA_character_
            )[1]
            multPairs <- x[grep("aligned concordantly >1 times$", x)]
            multPairs <- c(
                gsub("([0-9]*) .+ >1 times", "\\1", multPairs), NA_character_
            )[1]
            uniqDiscord <- x[grep("aligned discordantly 1 time$", x)]
            uniqDiscord <- c(
                gsub("([0-9]*) .+ 1 time", "\\1", uniqDiscord), NA_character_
            )[1]
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
        ## Add the Alignment Rate as a numeric
        out$Alignment_Rate <- alnRate
        as.data.frame(out)
    })

    ## Bind all result together
    df <- bind_rows(df)

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
        cols <- gsub(
            "([ACGTacgt]{2}/[ACGTacgt]{2})", "\\U\\1", cols, perl = TRUE
        )
        ## Form a list then tibble
        out <- as.list(x[,2])
        names(out) <- cols
        ## Using tibble here preserves column names in a nicer format
        as_tibble(out)
    })
    df <- dplyr::bind_rows(df)

    ## Reformat the time columns
    timeCols <- grepl("On$", names(df))
    df[timeCols] <- lapply(df[timeCols], parse_date_time, orders = "b! d! HMS")

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
        "Filename", starts_with("Total"),
        contains("Input"),
        contains("Mapped"), contains("Splice"),
        ends_with("On"), contains("Mapping"),
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
        which <- match.arg(which, c("metrics", "histogram", "all"))
    if (is.numeric(which)) stopifnot(which %in% c(1, 2))

    ## Collect the metrics from all files as a tibble
    metrics <- lapply(data, function(x){
        # ## Find the library name
        cmd <- stringr::str_subset(x, "(INPUT|input)")
        ip <- stringr::str_replace_all(
            cmd, ".+(input|INPUT)([ =]*)([^ ]+.bam).+", "\\3"
        )
        ip <- str_replace_all(ip, "[\\[\\]]", "")
        ## Find the header. The next two rows will be the colnames + data
        metHeader <- grep("METRICS CLASS\t.+.DuplicationMetrics", x)
        df <- .splitByTab(x[seq(metHeader + 1, by = 1, length.out = 2)])
        df$LIBRARY <- ip
        df
    })
    metrics <- dplyr::bind_rows(metrics)
    ## Now ensure the correct types
    metrics$PERCENT_DUPLICATION <- as.numeric(metrics$PERCENT_DUPLICATION)
    intCols <- setdiff(colnames(metrics), c("LIBRARY", "PERCENT_DUPLICATION"))
    metrics[intCols] <- lapply(metrics[intCols], as.integer)
    metrics <- as_tibble(metrics)

    ## Collect the histogram data from all files as a tibble
    ## This may not be present, so this is an optional step
    histData <- lapply(data, function(x){
        # ## Find the library name
        cmd <- stringr::str_subset(x, "(INPUT|input)")
        ip <- stringr::str_replace_all(
            cmd, ".+(input|INPUT)([ =]*)([^ ]+.bam).+", "\\3"
        )
        ip <- str_replace_all(ip, "[\\[\\]]", "")

        ## Find the header then remove up until that line
        histHeader <- grep("HISTOGRAM\tjava.lang.Double", x)
        df <- NULL
        if (length(histHeader)) {
            x <- x[-seq_len(histHeader)]
            ## Remove any blank lines (this is the last line)
            x <- x[!grepl("^$", x)]
            df <- .splitByTab(x)
            df$LIBRARY <- ip
            df <- dplyr::select(df, "LIBRARY", everything())
        }
        df
    })
    histData <- dplyr::bind_rows(histData)
    ## Ensure the correct types
    if (ncol(histData)) {
        histData <- lapply(
            histData,
            function(x) {
                any_na <- any(is.na(suppressWarnings(as.numeric(x))))
                if (!any_na) x <- as.numeric(x)
                x
            }
        )
        histData$BIN <- as.integer(histData$BIN)
        histData <- as_tibble(histData)
    }

    ## Setup the output, then format the column names
    out <- list(metrics = metrics, histogram = histData)
    out <- lapply(
        out,
        function(x){
            colnames(x) <- stringr::str_replace_all(colnames(x), "_", " ")
            colnames(x) <- stringr::str_to_title(colnames(x))
            x
        })
    if (which != "all") out <- out[[which]]
    out
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
    modNames <- c("sequences", "settings", "statistics", "distribution", "all")
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
            string = x[vals2Mods == 1][-1], pattern = ": ", n =  2
        )
        sequences <- as_tibble(
            split(sequences[,2], f = sequences[,1])
        )
        ## Now the settings
        settings <- str_split_fixed(
            string = x[vals2Mods == 2][-1], pattern = ": ", n =  2
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
            NA_integer_, as.integer(settings[["RNG seed"]])
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
            sequences = sequences, settings = settings,
            statistics = statistics, distribution = distribution
        )
        if (which != 'all') out <- out[which]
        out
    }

    arOut <- lapply(data, .parseArSingle, which = which)
    nMod <- unique(vapply(arOut, length, integer(1)))
    stopifnot(length(nMod) == 1)
    arOut <- lapply(
        seq_len(nMod),
        function(i) {
            dplyr::bind_rows(
                lapply(arOut, function(x) x[[i]]), .id = "Filename"
            )
        }
    )
    if (which == 'all') {
        names(arOut) <- modNames[seq_len(4)]
    } else {
        arOut <- arOut[[1]]
    }
    arOut

}

#' @title Parse data from cutadapt log files
#' @description Parse data from cutadapt log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param which which element of the log file to return.
#' Can be summary, adapter1, adapter2, adapter3 or overview, or any integer in
#' 1:5
#' @return tibble
#' @keywords internal
#' @importFrom tidyselect everything
#' @importFrom dplyr mutate_at mutate_if
#' @importFrom stringr str_replace_all str_remove_all str_extract
.parseCutadaptLogs <- function(data, which = 1){

    ## Possible module names vary depending on PE/SE reads and the adapter
    ## configuration. The only stable name is `summary`.
    ## For this reason module name checking will be handled inside the
    ## parsing function
    isMinimal <- vapply(data, length, integer(1)) == 2
    if (sum(isMinimal) > 0 & mean(isMinimal) != 1) stop(
        "A mix of minimal and full cutadapt reports was detected"
    )
    if (all(isMinimal)) which <- 1

    .parseCutAdaptSingle <- function(x, which){

        ## Check if it is a minimal report & return the df if TRUE
        if (length(x) == 2) {
            df <- .splitByTab(x)
            ## Values may sometimes be > 2^31 so convert to doubles
            df <- dplyr::mutate_if(df, colnames(df) != "status", as.numeric)
            return(df)
        }

        ## Otherwise parse the modules:
        ## Remove the blank elements
        x <- x[x != ""]
        ## Find where the modules start
        mods <- grepl("===", x)
        foundMods <- tolower(str_remove_all(x[mods],"[= ]"))

        ## Check the requested module
        if (is.numeric(which) & which > length(foundMods)) stop(
            "Invalid module selected. The maximum number is ",
            length(foundMods)
        )
        if (is.character(which) & which != "all") {
            if (sum(grepl(which, foundMods)) != 1) stop(
                "Unable to determine module. Valid modules are ",
                paste(foundMods, collapse = "/")
            )
        }

        ## Now split into the modules
        x <- split(x, f = cumsum(mods))
        names(x) <- c("header", foundMods)
        ## Grab the header as a list
        hdr <- x[["header"]]
        names(hdr) <- c("version", "params", "processing", "finished")
        ## Remove the first value from each now we have the intact header
        x <- lapply(x, function(x){x[-1]})

        ## Just do the modules manually
        out <- vector("list", length(x) - 1)
        names(out) <- names(x)[-1]

        ## Start with the summary (i.e. minimal). Values then names
        vals <- str_replace_all(x$summary, ".+ +([0-9,]+).+", "\\1")
        vals <- str_remove_all(vals, ",")
        vals <- as.numeric(vals)
        nm <- str_trim(
            str_replace_all(x$summary, "(.+):.+", "\\1")
        )
        ## If this is paired end, some of these will be nested and need to be
        ## handled correctly. The first two will be processed bp, whilst
        ## the next two will be written bp. Using the above code, they should
        ## just appear as Read 1 and Read 2 for each set.
        dups <- grepl("^Read [12]$", nm)
        nm[dups] <- paste(
            nm[dups], "basepairs",
            rep(c("processed", "written"), each = sum(dups)/2)
        )
        nm <- tolower(nm)
        nm <- str_remove_all(nm, "[\\(\\)]")
        nm <- str_replace_all(nm, " ", "_")
        names(vals) <- nm
        out[["summary"]] <- as_tibble(as.list(vals))
        out[["summary"]][["header"]] <- list(hdr)
        out[["summary"]] <- dplyr::select(
            out[["summary"]], "header", everything()
        )

        ## Now the adapters
        out[grep("adapter", names(out))] <- lapply(
            x[grepl("adapter", names(x))],
            function(a){
                ## Avoid R CMD errors by declaring dummy variables
                A <- C <- G <- `T` <- c()
                main <- tibble(
                    sequence = str_extract(a[grepl("Sequence", a)], "[ACGTN]+"),
                    type = str_replace_all(
                        a[grepl("Type", a)], ".+Type: (.+); Length.+", "\\1"
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
                    ) / 100,
                    `none/other` = round(1 - (A + C + G + `T`), 2)
                )
                overview <- .splitByTab(
                    a[seq(which(grepl("^Overview", a)) + 1, length(a))]
                )
                ## Needs to be recoded using `across`
                overview <- dplyr::mutate_if(
                    overview,
                    vapply(
                        overview,
                        \(x) suppressWarnings(all(!is.na(as.numeric(x)))),
                        logical(1)
                    ),
                    as.numeric
                )
                main[["overview"]] <- list(as_tibble(overview))
                main
            }
        )

        ## Return the required module
        if (which != "all") out <- out[which]
        out
    }

    caOut <- lapply(data, .parseCutAdaptSingle, which = which)
    if (all(isMinimal)) {
        caOut <- bind_rows(caOut, .id = "Filename")
        return(caOut)
    }

    nMod <- unique(vapply(caOut, length, integer(1)))
    stopifnot(length(nMod) == 1)
    modNames <- as.character(unlist(lapply(caOut, names)))
    caOut <- lapply(
        seq_len(nMod),
        \(i) dplyr::bind_rows(
            lapply(caOut, function(x) x[[i]]), .id = "Filename"
        )
    )
    if (which == 'all') {
        names(caOut) <- modNames
    } else {
        caOut <- caOut[[1]]
    }
    caOut

}

#' @title Parse data from featureCounts summary files
#' @description Parse data from featureCounts summary files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return tibble
#' @keywords internal
#' @importFrom forcats fct_inorder
#' @importFrom tidyselect ends_with
.parseFeatureCountsLogs <- function(data, ...){

    out <- lapply(data, function(x){
        x <- .splitByTab(x)
        x <- tidyr::pivot_longer(
            data = x,
            cols = ends_with("bam"),
            names_to = "Sample",
            values_to = "Total"
        )
        x$Sample <- x$Sample
        x$Total <- as.integer(x$Total)
        x$Status <- fct_inorder(x$Status)
        tidyr::pivot_wider(
            data = x,
            id_cols = "Sample",
            names_from = "Status",
            values_from = "Total"
        )
    })

    out <- lapply(names(data), function(x){
        out[[x]]$Filename <- x
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
        ## get numbers of each field from lines
        single <- grep(
            "Complete and single-copy BUSCOs \\(S\\)", x = x, value = TRUE
        )
        single <- as.integer(gsub(".*\t(.+)\tComp.*", "\\1", single))

        duplicated <- grep(
            "Complete and duplicated BUSCOs \\(D\\)", x = x, value = TRUE
        )
        duplicated <- as.integer(gsub(".*\t(.+)\tComp.*", "\\1", duplicated))

        fragmented <- grep("Fragmented BUSCOs \\(F\\)", x = x, value = TRUE)
        fragmented <- as.integer(gsub(".*\t(.+)\tFrag.*", "\\1", fragmented))

        missing <- grep("Missing BUSCOs \\(M\\)", x = x, value = TRUE)
        missing <- as.integer(gsub(".*\t(.+)\tMiss.*", "\\1", missing))

        total <- sum(single, duplicated, fragmented, missing)

        name <- grep(
            "Summarized benchmarking in BUSCO notation for file ",
            x = x, value = TRUE
        )
        name <- gsub(".*file ", "\\1", name)
        df <- tibble(
            name = name,
            completeSingleCopy = single,
            completeDuplicated = duplicated,
            fragmented = fragmented, missing = missing
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

            Input_Read_Pairs <- gsub(
                "Input Read Pairs: ([0-9]+).+", "\\1", valLine
            )
            Input_Read_Pairs <- as.integer(Input_Read_Pairs)

            Both_Surviving <- gsub(
                ".+Both Surviving: ([0-9]+).+", "\\1", valLine
            )
            Both_Surviving <- as.integer(Both_Surviving)

            Forward_Only_Surviving <- gsub(
                ".+Forward Only Surviving: ([0-9]+).+", "\\1", valLine
            )
            Forward_Only_Surviving <- as.integer(Forward_Only_Surviving)

            Reverse_Only_Surviving <- gsub(
                ".+Reverse Only Surviving: ([0-9]+).+", "\\1", valLine
            )
            Reverse_Only_Surviving <- as.integer(Reverse_Only_Surviving)
        }
        Dropped <- gsub(".+Dropped: ([0-9]+).+", "\\1", valLine)
        Dropped <- as.integer(Dropped)

        tibble(
            Type = readType,
            Input_Reads, Input_Read_Pairs,
            Surviving, Both_Surviving,
            Forward_Only_Surviving, Reverse_Only_Surviving, Dropped,
            Illumina_Clip, Sliding_Window, Max_Info,
            Leading, Trailing,
            Crop, Head_Crop,
            Min_Len, Avg_Qual, Quality_Encoding
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
    out <- pivot_wider(data = out, names_from = "key", values_from = "value")

    ## Return the final output
    dplyr::select(
        out,
        "Filename", "Type", starts_with("Input"), contains("Surviving"),
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
            N50 = N50, N75 = N75,
            L50 = L50, L75 = L75
        )
    })
    df <- dplyr::bind_rows(df)

    dplyr::bind_cols(tibble(fileNames = names(data)), df)

}

#' @title Parse data from samtools flagstat files
#' @description Parse data from samtools flagstat files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseFlagstatLogs <- function(data, ...){

    .parseFlagstatSingle <- function(data){
        df <- gsub("([0-9]+) \\+ ([0-9]+) (.+)", "\\1\t\\2\t\\3", data)
        df <- .splitByTab(df, FALSE)
        colnames(df) <- c("QC-passed", "QC-failed", "flag")
        df[["flag"]] <- gsub("([a-z ]+) \\(.+", "\\1", df[["flag"]])
        df[["flag"]][nrow(df)] <- paste(df[["flag"]][nrow(df)], "(mapQ>=5)")
        df[["QC-passed"]] <- as.integer(df[["QC-passed"]])
        df[["QC-failed"]] <- as.integer(df[["QC-failed"]])
        df
    }

    out <- lapply(data, .parseFlagstatSingle)
    out <- lapply(names(out), function(x){
        out[[x]]$Filename <- x
        out[[x]]
    })
    out <- bind_rows(out)
    dplyr::select(out, "Filename", everything())


}

#' @title Parse data from macs2 callpeak log files
#' @description Parse data from macs2 callpeak log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @importFrom lubridate as_datetime
#' @importFrom stringr str_subset str_remove_all str_split_fixed str_extract
#' @importFrom stringr str_replace_all str_split str_trim str_extract_all
#' @importFrom tidyselect everything ends_with contains any_of
#' @keywords internal
.parseMacs2CallpeakLogs <- function(data, ...){

    .parseSingleCallPeak <- function(data, ...){

        ## Remove any leading/trailing whitespace
        data <- str_trim(data)

        ## The date the command was run
        dt <- gsub("INFO.+@ (.+):", "\\1", data[[1]])
        dt <- as_datetime(
            dt, format = "%a, %d %b %Y %H:%M:%S", tz = Sys.timezone()
        )

        ## Fragment Length
        fl <- data[grepl("INFO .+ fragment (length|size) is", data) ]
        fl <- gsub(".+fragment (length|size) is.* ([0-9\\.]+) bp.*", "\\2", fl)
        fl <- as.numeric(fl)

        ## Tag length
        tl <- data[grepl("INFO .+ tag size = ", data) ]
        tl <- gsub(".+tag size = ([0-9\\.]+)$", "\\1", tl)
        tl <- as.numeric(tl)

        ## The arguments in the main header
        ind <- seq(
            which(grepl("ARGUMENTS LIST", data)),
            min(which(data == "")) - 1
        )
        args <- data[ind]
        args_out <- args[grepl(" = ", args)]
        args_out <- str_remove_all(args_out, "# ")
        args_out <- str_split_fixed(args_out, " = ", 2)
        args_out[, 1] <- gsub("[ -]", "_", args_out[, 1])
        args_out <- structure(args_out[, 2], names = args_out[, 1])
        args_out <- as.list(args_out)
        nums <- suppressWarnings(vapply(args_out, as.numeric, numeric(1)))
        args_out[!is.na(nums)] <- lapply(args_out[!is.na(nums)], as.numeric)
        args_out$date <- dt

        ## Add any additional arguments
        cmd <- data[grepl("Command line: callpeak", data)]
        if (!grepl("--max-gap", cmd)) args_out$max_gap <- tl
        if (!grepl("--min-length", cmd)) args_out$min_length <- fl
        args_out$keep_dup <- grepl("--keep-dup", cmd)
        args_out$nomodel <- grepl("nomodel", cmd)

        args_out$scale_to <- str_subset(args, "dataset will be scaled towards")
        args_out$scale_to <- str_extract(args_out$scale_to , "(large|small)")

        lambda <- str_subset(args, "regional lambda")
        lambda <- str_extract_all(lambda, "[0-9]+ bps")
        lambda <- lapply(lambda, str_remove_all, " bps")
        lambda <- rep(unlist(lambda), 2)[seq_len(2)]
        lambda <- paste0("[", paste(lambda, collapse = ", "), "]")
        args_out$local <- lambda

        args_out$broad <- grepl("--broad ", cmd)
        args_out$paired_end <- str_subset(args, "Paired-End mode")
        args_out$paired_end <- !grepl("off", args_out$paired_end)

        if ("ChIP_seq_file" %in% names(args_out)) {
            l <- args_out$ChIP_seq_file
            l <- str_remove_all(l, "[\\[\\]\\']")
            l <- str_split(l, pattern = ", ")[[1]]
            args_out$ChIP_seq_file <- list(l)
            rm(l)
        }

        if ("control_file" %in% names(args_out)) {
            l <- args_out$control_file
            l <- str_remove_all(l, "[\\[\\]\\']")
            l <- str_split(l, pattern = ", ")[[1]]
            args_out$control_file <- list(l)
            rm(l)
        }

        ## Read numbers
        n_reads <- str_subset(data, "(read|fragment)s have been read")
        n_reads <- str_replace_all(
            n_reads, "INFO.+ ([0-9]+) (read|fragment)s have been.+", "\\1"
        )
        n_reads <- as.numeric(n_reads)
        args_out$n_reads <- list(n_reads)

        ## Tag numbers
        n_tags_treatment <- str_subset(data, "total tags in treatment")
        n_tags_treatment <- str_replace_all(
            n_tags_treatment , ".+ ([0-9]+)$", "\\1"
        )
        args_out$n_tags_treatment <- as.numeric(n_tags_treatment)
        n_tags_control <- str_subset(data, "total tags in control")
        n_tags_control <- str_replace_all(
            n_tags_control, ".+ ([0-9]+)$", "\\1"
        )
        args_out$n_tags_control <- as.numeric(n_tags_control)

        ## Peaks
        paired_peaks <- str_subset(data, "number of paired peaks")
        paired_peaks <- str_replace_all(paired_peaks, ".+ ([0-9]+)$", "\\1")
        args_out$paired_peaks <- as.numeric(paired_peaks)

        ## Fragment length & tag length
        args_out$tag_length <- tl
        args_out$fragment_length <- fl
        alt <- str_subset(data, "alternative fragment length")
        alt <- str_replace_all(alt, ".+ ([0-9,]+) bps.*", "\\1")
        alt <- str_split(alt, ",")
        alt <- lapply(alt, as.numeric)
        args_out$alt_fragment_length <- alt

        ## Output files
        r_script <- str_subset(data, "R script")
        r_script <- str_replace_all(
            r_script, ".+R script for model : (.+)", "\\1"
        )
        xls <- str_subset(data, "Write output xls file")
        xls <- str_replace_all(xls, ".+Write output xls file.+ (.+)", "\\1")

        peakFiles <- str_subset(data, " Write .* in .+Peak format file")
        peakFiles <- str_remove_all(peakFiles, ".+format file\\.\\.\\. ")

        summits_bed <- str_subset(data, " Write summits bed file")
        summits_bed <- str_replace_all(
            summits_bed, ".+ Write summits bed file.+ (.+)", "\\1"
        )
        args_out$outputs <- list(
            c(
                r_script, xls, peakFiles, summits_bed
            )
        )

        ## To avoid R CMD check issues
        name <- min_length <- c()
        ## Suggested by LaureTomas
        i <- vapply(args_out, length, integer(1)) > 0
        args_out <- args_out[i]
        args_out <- as_tibble(args_out)
        dplyr::select(
            args_out,
            name, date, any_of("paired_peaks"), min_length, contains("tags"), n_reads,
            ends_with("length"), everything()
        )
    }

    ## Run using lapply, then bind_rows and add the Filenames
    df <- lapply(data, .parseSingleCallPeak)
    df <- dplyr::bind_rows(df)
    df$Filename <- names(data)
    dplyr::select(df, "Filename", everything())


}

#' @title Parse data from umi_tools dedup files
#' @description Parse data from umi_tools dedup files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @importFrom stats setNames
#' @keywords internal
.parseUmitoolsDedupLogs <- function(data, ...){

    .parseSingle <- function(data){
        timing <- data[[length(data)]]
        liveLines <- grepl("(Written out|Parsed|Searching)", data)
        commentLines <- grepl("^#", data)
        warningLines <- grepl("WARNING", data)
        cmd <- grepl("command", data)
        data <- data[!(liveLines | commentLines | warningLines | cmd)]
        data <- gsub("(.+INFO )", "", data)
        data <- gsub("^Reads: ", "", data)
        data <- strsplit(data, ",")
        data <- unlist(lapply(data, strsplit, ":"))
        data <- stringr::str_trim(data)
        mates <- data[grepl("mates", data)]
        data <- data[!grepl("mates", data)]
        umis <- data[grepl("umis", data)]
        data <- data[!grepl("umis", data)]
        mat <- matrix(data, ncol = 2, byrow = TRUE)
        named_vec <- setNames(mat[,2], gsub("[ \\.]+", "_", mat[,1]))
        named_list <- as.list(named_vec)
        ## The next two fields were introduced with v0.23
        ## Should just duplicate the number of reads
        total_umis <- gsub("total_umis ", "", umis[grepl("total_umis", umis)])
        if (length(total_umis)) named_list$total_umis <- total_umis
        ## The total # unique UMIs
        unique_umis <- gsub("#umis ", "", umis[grepl("#umis ", umis)])
        if (length(unique_umis)) named_list$unique_umis <- unique_umis
        named_list$Mates_never_found <- gsub(" mates never found", "", mates)
        named_list$start <- gsub(".+at (.+) --.+--.+", "\\1", timing)
        named_list$duration <- gsub(".+finished in (.+) at.+", "\\1", timing)
        as_tibble(named_list)
    }

    df_list <- lapply(data, .parseSingle)
    df <- dplyr::bind_rows(df_list, .id = "Filename")
    num_cols <- grepl("Mean", colnames(df))
    df[num_cols] <- lapply(df[num_cols], as.numeric)
    int_cols <- grepl("([Rr]ead|Total|Max|Mates)", colnames(df))
    df[int_cols] <- lapply(df[int_cols], as.integer)
    df$start <- lubridate::as_datetime(
        df$start, format = "%a %b %d %H:%M:%S %Y", tz = Sys.timezone()
    )
    df$duration <- lubridate::as.duration(df$duration)
    df

}

#' @title Parse data from rnaseqc metrics files
#' @description Parse data from rnaseqc metrics files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @importFrom stats setNames
#' @keywords internal
.parseRnaseqcMetricsLogs <- function(data, ...){

    df_list <- lapply(data, .splitByTab, firstRowToNames = FALSE)
    tbl_list <- lapply(df_list, \(x) setNames(as.list(x[[2]]), x[[1]]))
    tbl_list <- lapply(tbl_list, tibble::as_tibble) # Respects column names
    tbl <- dplyr::bind_rows(tbl_list, .id = "Filename")
    nm <- names(tbl)
    i <- setdiff(nm, c("Sample", "Filename"))
    tbl[i] <- lapply(tbl[i], as.numeric)
    tbl

}

#' @title Parse data from salmon lib_format_counts.json files
#' @description Parse data from salmon lib_format_counts.json files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseSalmonLibFormatCountsLogs <- function(data, ...){

    js_list <- lapply(data, jsonlite::fromJSON)
    tbl_list <- lapply(js_list, tibble::as_tibble)
    tbl <- dplyr::bind_rows(tbl_list, .id = "Filename")
    tbl$read_files <- strsplit(
        gsub("(^\\[ |\\]$)", "", tbl$read_files), ", "
    )
    tbl
}

#' @title Parse data from salmon meta_info.json files
#' @description Parse data from salmon meta_info.json files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @param ... Not used
#' @return data.frame
#' @keywords internal
.parseSalmonMetaInfoLogs <- function(data, ...){

    js_list <- lapply(data, jsonlite::fromJSON)
    ## Remove any missing columns & tidy up other problematic columns
    js_list <- lapply(
        js_list,
        \(x) {
            ln <- vapply(x, length, integer(1))
            x[ln > 1] <- lapply(x[ln > 1], list)
            x[ln > 0]
        }
    )
    tbl_list <- lapply(js_list, tibble::as_tibble)
    tbl <- dplyr::bind_rows(tbl_list, .id = "Filename")
    time_cols <- grepl("time$", names(tbl))
    tbl[time_cols] <- lapply(
        tbl[time_cols], lubridate::as_datetime,
        format = "%a %b %d %H:%M:%S %Y", tz = Sys.timezone()
    )
    tbl
}
