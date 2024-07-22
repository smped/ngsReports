#' @importFrom methods as
#' @importFrom tidyselect one_of contains ends_with
setAs(".FastqcFile", "FastqcData", function(from){

    ## Import each line as an element in a character vector then detect the
    ## modules and return as a list containing each module as an element
    fqcLines <- .getFqcLines(from)

    ## The FastQC version number
    ## May not be present so is only returned if present
    vers <- NA_character_ # Default to missing
    nm <- tolower(names(fqcLines))
    if (grepl("fastqc\\t", nm[1])) {
        vers <- gsub("fastqc\\t(.+)", "\\1", nm[1])
        ## Remove the version module for easier module identification
        fqcLines <- fqcLines[-1]
    }

    ## Initialise an empty list with an empty element for each module present
    ## in the data, and the set of required modules
    out <- .initialiseOutput(fqcLines)

    ## Now wrangle the data for each module and populate the output
    ## Missing or empty modules will return an empty dataframe
    out[["Basic_Statistics"]] <- .getBasicStatistics(fqcLines)
    out[["Per_base_sequence_quality"]] <- .getPerBaseSeqQuals(fqcLines)
    out[["Per_tile_sequence_quality"]] <- .getPerTileSeqQuals(fqcLines)
    out[["Per_sequence_quality_scores"]] <- .getPerSeqQualScores(fqcLines)
    out[["Per_base_sequence_content"]] <- .getPerBaseSeqCont(fqcLines)
    out[["Per_sequence_GC_content"]] <- .getPerSeqGcCont(fqcLines)
    out[["Per_base_N_content"]] <-  .getPerBaseNCont(fqcLines)
    out[["Sequence_Length_Distribution"]] <- .getSeqLengthDist(fqcLines)
    out[["Overrepresented_sequences"]] <- .getOverrepSeq(fqcLines)
    out[["Adapter_Content"]] <- .getAdapterCont(fqcLines)
    out[["Kmer_Content"]] <- .getKmerCont(fqcLines)

    ## Get the Sequence Duplication Levels. This is potentially a two element
    ## module within the parsed data, depending on the version of FastQC
    seqDupLevels <- .getSeqDuplicationLevels(fqcLines)
    dupMods <- names(seqDupLevels)
    out[dupMods] <- seqDupLevels[dupMods]

    ## Get the summary information
    Summary <- getSummary(from)

    args <- list(
        Class = "FastqcData",
        path = path(from),
        version = vers,
        Summary = Summary
    )

    do.call("new", c(args, out))

})
#' @importFrom methods as
setAs("list", "FastqcDataList", function(from){
    isFqcData <- vapply(from, is, logical(1), class2 = "FastqcData")
    if (!all(isFqcData)) stop(
        "All elements of the supplied list must be FastqcData objects"
    )
    new("FastqcDataList", from)
})

#' @importFrom jsonlite validate fromJSON
setAs(".FastpFile", "FastpData", function(from){

    ## Import each line as an element in a character vector then detect the
    ## modules and return as a list containing each module as an element
    path <- path(from)
    lines <- suppressWarnings(readLines(path))
    stopifnot(validate(lines))
    data <- fromJSON(lines)
    exp_mods <- c(
        "summary", "filtering_result", "duplication", "adapter_cutting",
        "command"
    )
    stopifnot(all(exp_mods %in% names(data)))

    ## The fastp version number
    ## May not be present so is only returned if present
    vers <- NA_character_ # Default to missing
    if ("fastp_version" %in% names(data[["summary"]]))
        vers <- data[["summary"]][["fastp_version"]]

    ## Initialise an empty list too parse the data into
    ## Modules should be as defined in the man page for FastpData
    isPaired <- grepl("(-I|--in2)", data$command)
    out <- list()
    out[["Summary"]] <- .getFpSummary(data)
    out[["Adapters"]] <- .getFpAdapters(data)
    out[["Duplication"]] <- .getFpDuplication(data)
    out[["Insert_size"]] <- .getFpInsert(data)
    out[["Before_filtering"]] <- .getFpBeforeAfter(data, "before", isPaired)
    out[["After_filtering"]] <- .getFpBeforeAfter(data, "after", isPaired)
    out[["paired"]] <- isPaired
    out[["command"]] <- data[["command"]]
    out[["version"]] <- vers
    out[["path"]] <- path

    args <- list(Class = "FastpData")
    do.call("new", c(args, out))

})

.getFpSummary <- function(data){

    n_reads <- data$summary$before_filtering$total_reads

    ## The Before/After Filtering Table
    Before_filtering <- as_tibble(data$summary$before_filtering)
    After_filtering <- as_tibble(data$summary$after_filtering)

    ## The filtering Result
    results <- unlist(data$filtering_result)
    rates <- results  / n_reads
    filtering_result <- tibble(
        result = names(results), total = results,rate = rates
    )

    list(
        Before_filtering = Before_filtering, After_filtering = After_filtering,
        Filtering_result = filtering_result
    )
}

#' @importFrom stringr str_count
.getFpAdapters <- function(data){
    n_reads <- data$summary$before_filtering$total_reads
    n_bases <- data$summary$before_filtering$total_bases
    ad <- data$adapter_cutting
    ## Summary table
    tbl <- as_tibble(ad[!grepl("adapter_counts", names(ad))])
    tbl$adapter_trimmed_reads_rate <- tbl$adapter_trimmed_reads / n_reads
    tbl$adapter_trimmed_bases_rate <- tbl$adapter_trimmed_bases / n_bases

    ## R1
    r1 <- unlist(ad$read1_adapter_counts)
    tbl_r1 <- tibble(Sequence = names(r1), Occurences = as.integer(r1))
    tbl_r1$Occurence_rate <- tbl_r1$Occurences / n_reads
    tbl_r1$adapter_length <- str_count(tbl_r1$Sequence, "[ACGTN]")
    tbl_r1$adapter_length[tbl_r1$Sequence == "others"] <- NA_integer_

    ## R2
    r2 <- unlist(ad$read2_adapter_counts)
    tbl_r2 <- tibble(Sequence = names(r2), Occurences = as.integer(r2))
    tbl_r2$Occurence_rate <- tbl_r2$Occurences / n_reads
    tbl_r2$adapter_length <- str_count(tbl_r2$Sequence, "[ACGTN]")
    tbl_r2$adapter_length[tbl_r2$Sequence == "others"] <- NA_integer_
    tbl_r2

    tbl$read1_adapter_count <- list(tbl_r1)
    tbl$read2_adapter_count <- list(tbl_r2)
    tbl
}

.getFpDuplication <- function(data){
    tbl <- tibble(rate = data$duplication$rate)
    hist <- tibble(
        histogram = data$duplication$histogram,
        mean_gc = data$duplication$mean_g,
    )
    hist$duplication_rate <- numeric(nrow(hist))
    hist$duplication_level <- seq_len(nrow(hist))
    if (nrow(hist))
        hist$duplication_rate <- hist$histogram / sum(hist$histogram)
    tbl$histogram <- list(hist)
    tbl
}

#' @importFrom tidyr nest
.getFpInsert <- function(data){
    insert_size <- histogram <- freq <- NULL
    n_sampled <- c(data$insert_size$unknown, data$insert_size$histogram)
    tbl <- as_tibble(data$insert_size)
    tbl$unknown_rate <- tbl$unknown / sum(n_sampled)
    tbl$freq <- tbl$histogram / sum(n_sampled)
    tbl$insert_size <- seq_along(tbl$histogram)
    nest(tbl, histogram = c(insert_size, histogram, freq))
}

.getFpBeforeAfter <- function(data, type = c("before", "after"), paired){

    type <- match.arg(type)
    count <- NULL

    ## Setup the kmer values
    atcg <- c("A", "T", "C", "G")
    col_levels <- apply(
        expand.grid(atcg, atcg)[,c(2,1)], 1, paste, collapse = ""
    )
    row_levels <- apply(
        expand.grid(atcg, atcg, atcg)[,c(3, 2,1)], 1, paste, collapse = ""
    )
    cols <- c(
        "total_reads", "total_bases", "q20_bases", "q30_bases", "total_cycles"
    )

    ## Read1 should always be present
    mod <- paste(c("read1_", type, "_filtering"), collapse = "")
    stopifnot(mod %in% names(data))
    tbl <- as_tibble(data[[mod]][cols])

    quality_curves <- as_tibble(data[[mod]]$quality_curves)
    quality_curves$position <- seq_len(nrow(quality_curves))
    tbl$quality_curves <- list(quality_curves)

    content_curves <- as_tibble(data[[mod]]$content_curves)
    content_curves$position <- seq_len(nrow(content_curves))
    tbl$content_curves <- list(content_curves)

    kmer_count <- unlist(data[[mod]]$kmer_count)
    kmer_tbl <- tibble(
        kmer = names(kmer_count), count = as.integer(kmer_count)
    )
    kmer_tbl$times_mean <- kmer_tbl$count / mean(kmer_tbl$count)
    kmer_tbl$prefix <- gsub("^([ATCG]{3}).+", "\\1", kmer_tbl$kmer)
    kmer_tbl$prefix <- factor(kmer_tbl$prefix, levels = row_levels)
    kmer_tbl$suffix <- gsub("^([ATCG]{3})([ATCG]{2})", "\\2", kmer_tbl$kmer)
    kmer_tbl$suffix <- factor(kmer_tbl$suffix, levels = col_levels)
    tbl$kmer_count <- list(kmer_tbl)

    overrep <- data[[mod]]$overrepresented_sequences
    if (length(overrep) > 0) {
        tbl$overrepresented_sequences <- list(
            tibble(
                sequence = names(overrep),
                count = as.integer(unlist(overrep)),
                ## I don't know why we multiply by 2, but it gives the same
                ## value as the html report
                rate = nchar(sequence) * count * 2 / tbl$total_bases
            )
        )
    }

    out <- list(read1 = tbl)

    ## Check for the read2 module
    if (paired) {
        mod <- paste(c("read2_", type, "_filtering"), collapse = "")
        stopifnot(mod %in% names(data))
        tbl <- as_tibble(data[[mod]][cols])

        quality_curves <- as_tibble(data[[mod]]$quality_curves)
        quality_curves$position <- seq_len(nrow(quality_curves))
        tbl$quality_curves <- list(quality_curves)

        content_curves <- as_tibble(data[[mod]]$content_curves)
        content_curves$position <- seq_len(nrow(content_curves))
        tbl$content_curves <- list(content_curves)

        kmer_count <- unlist(data[[mod]]$kmer_count)
        kmer_tbl <- tibble(
            kmer = names(kmer_count), count = as.integer(kmer_count)
        )
        kmer_tbl$times_mean <- kmer_tbl$count / mean(kmer_tbl$count)
        kmer_tbl$prefix <- gsub("^([ATCG]{3}).+", "\\1", kmer_tbl$kmer)
        kmer_tbl$prefix <- factor(kmer_tbl$prefix, levels = row_levels)
        kmer_tbl$suffix <- gsub("^([ATCG]{3})([ATCG]{2})", "\\2", kmer_tbl$kmer)
        kmer_tbl$suffix <- factor(kmer_tbl$suffix, levels = col_levels)
        tbl$kmer_count <- list(kmer_tbl)

        overrep <- data[[mod]]$overrepresented_sequences
        if (length(overrep) > 0) {
            tbl$overrepresented_sequences <- list(
                tibble(
                    sequence = names(overrep),
                    count = as.integer(unlist(overrep)),
                    ## I don't know why we multiply by 2, but it gives the same
                    ## value as the html report
                    rate = nchar(sequence) * count * 2 / tbl$total_bases
                )
            )
        }
        out[["read2"]] <- tbl
    }
    out
}


## This helper checks for compressed or extracted FastQC reports then
## imports the contents of fastqc_data.txt as a character vector.
## Comments (#) and lines denoting the end of a module are then removed
.getFqcLines <- function(x){

    ## Setup the path to the file & use tryCatch to ensure it exists
    ## when testing to see if it is compressed
    path <- path(x)
    stopifnot(file.exists(path))
    comp <- isCompressed(path, type = "zip")

    if (comp) {

        ## Get the list of files within the zip archive
        allFiles <- unzip(path, list = TRUE)$Name
        ## Check 'fastqc_data.txt' exists within the file
        subdir <- gsub(".zip$", "", basename(allFiles[1]))
        fl <- file.path(subdir, "fastqc_data.txt")
        stopifnot(fl %in% allFiles)

        ## # Open the connection & read the 12 lines
        uz <- unz(path, fl)
        fqcLines <- readLines(uz)
        close(uz)

    }
    else{

        ## The existence of this file will have been checked at
        ## instantiation of the .FastqcFile. Check in case it has been
        ## deleted post-instantiation
        fl <- file.path(path, "fastqc_data.txt")
        if (!file.exists(fl)) stop("'fastqc_data.txt' could not be found.")
        fqcLines <- readLines(fl)

    }

    ## Modified from Repitools::readFastQC
    ## Remove any '#' symbols
    fqcLines <- gsub("#", "", fqcLines)
    ## Remove any lines which specify '>>END_MODULE'
    fqcLines <- fqcLines[!grepl(">>END_MODULE", fqcLines)]

    ## Split the data based on the '>>' pattern, which indicates the
    ## beginning of a new module
    fqcLines <- split(fqcLines, cumsum(grepl("^>>", fqcLines)))

    ## Assign the module names
    names(fqcLines) <- vapply(fqcLines, function(x){
        ## Get the first line before the tab separator
        nm <- gsub(">>(.+)\\t.+", "\\1", x[[1]])
        ## Replace the space with underscore
        gsub(" ", "_", nm)
    }, character(1))
    ## Remove the module name (i.e. the 1st value) from each module's data
    fqcLines <- lapply(fqcLines, function(x){x[-1]})

    fqcLines

}

## This helper checks for the presence of required modules then initialises
## the empty list which will hold the module data
.initialiseOutput <- function(x){

    ## Get the module names from the lines that are present in the data
    modules <- names(x)

    ## Define the standard modules in a FastQC report. These are required
    reqModules <- c(
        "Basic_Statistics",
        "Per_base_sequence_quality",
        "Per_tile_sequence_quality",
        "Per_sequence_quality_scores",
        "Per_base_sequence_content",
        "Per_sequence_GC_content",
        "Per_base_N_content",
        "Sequence_Length_Distribution",
        "Sequence_Duplication_Levels",
        "Overrepresented_sequences",
        "Adapter_Content",
        "Kmer_Content"
    )

    ## Check that at least one of the standard modules is present
    ## It's quite plausible for a report to be missing multiple modules
    if (!any(modules %in% reqModules))
        stop("None of the standard modules were found in the data.")

    ## Define the output to have the same structure as fastqcData,
    ## with an additional slot to account for the Sequence Duplication
    ## Levels output. Initialise an empty list based on the standard
    ## modules after checking for the Total_Deduplicated Percentage value.
    allMods <- c(reqModules, modules, "Total_Deduplicated_Percentage")
    allMods <- unique(allMods)
    out <- vector("list", length(allMods))
    names(out) <- allMods

    ## Now return the initialised output
    out
}

# Define a series of quick functions for wrangling the data
# after splitting the input from readLines()
# These are all tested using testthat to ensure they're working correctly
.getBasicStatistics <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Basic_Statistics"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    ## Form a data frame and check th structure
    x <- .splitByTab(fqcLines[[mod]])
    stopifnot(setequal(names(x), c("Measure", "Value")))
    vals <- x[["Value"]]
    names(vals) <- gsub(" ", "_", x[["Measure"]])

    ## Check for the required values
    reqVals <- c(
        "Filename",
        "File_type",
        "Encoding",
        "Total_Sequences",
        "Sequences_flagged_as_poor_quality",
        "Sequence_length",
        "%GC"
    )
    stopifnot(reqVals %in% names(vals))

    ## Setup the tibble with the correct value types
    df <- as_tibble(as.list(vals))
    df$Shortest_sequence <- gsub("(.*)-.*", "\\1", df$Sequence_length)
    df$Longest_sequence <- gsub(".*-(.*)", "\\1", df$Sequence_length)
    intVals <- c(
        "Shortest_sequence",
        "Longest_sequence",
        "Total_Sequences",
        "Sequences_flagged_as_poor_quality"
    )
    df[intVals] <- lapply(df[intVals], as.integer)
    dplyr::select(
        df,
        "Filename",
        "Total_Sequences",
        contains("quality"),
        ends_with("sequence"),
        one_of("%GC", "File_type", "Encoding")
    )

}

.getPerBaseSeqQuals <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Per_base_sequence_quality"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])
    names(df) <- gsub(" ", "_", names(df))

    ## Check for the required values
    reqVals <- c(
        "Base",
        "Mean",
        "Median",
        "Lower_Quartile",
        "Upper_Quartile",
        "10th_Percentile",
        "90th_Percentile"
    )
    stopifnot(reqVals %in% names(df))

    ## Change to numeric where appropriate
    df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
    as_tibble(df)

}

.getPerTileSeqQuals <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Per_tile_sequence_quality"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])

    ## Check for the required values
    reqVals <- c("Tile", "Base", "Mean")
    stopifnot(reqVals %in% names(df))

    df[["Mean"]] <- as.numeric(df[["Mean"]])
    as_tibble(df)

}

.getPerSeqQualScores <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Per_sequence_quality_scores"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])

    ## Check for the required values
    reqVals <- c("Quality", "Count")
    stopifnot(reqVals %in% names(df))

    df <- lapply(df, as.integer)
    as_tibble(df)
}

.getPerBaseSeqCont <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Per_base_sequence_content"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])

    ## Check for the required values
    reqVals <- c("Base", "G", "A", "T", "C")
    stopifnot(reqVals %in% names(df))

    ## Set all as numeric except the 'Base'column
    df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
    as_tibble(df)

}

.getPerSeqGcCont <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Per_sequence_GC_content"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])
    names(df) <- gsub(" ", "_", names(df))

    ## Check for the required values
    reqVals <- c("GC_Content", "Count")
    stopifnot(reqVals %in% names(df))

    ## Convert to integers
    df[["GC_Content"]] <- as.integer(df[["GC_Content"]])
    df[["Count"]] <- as.numeric(df[["Count"]])
    as_tibble(df)
}

.getPerBaseNCont <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Per_base_N_content"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])

    ## Check for the required values
    reqVals <- c("Base", "N-Count")
    stopifnot(reqVals %in% names(df))

    df[["N-Count"]] <- as.integer(df[["N-Count"]])
    as_tibble(df)

}

.getSeqLengthDist <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Sequence_Length_Distribution"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])

    ## Check for the required values
    reqVals <- c("Length", "Count")
    stopifnot(reqVals %in% names(df))

    ## If there is a range, separate into Lower & Upper values
    df$Lower <- as.integer(gsub("(.*)-.*", "\\1", df$Length))
    df$Upper <- as.integer(gsub(".*-(.*)", "\\1", df$Length))
    df$Count <- as.integer(df$Count)
    df <- as_tibble(df)

    df[c("Length", "Lower", "Upper", "Count")]

}

.getSeqDuplicationLevels <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Sequence_Duplication_Levels"
    if (!mod %in% names(fqcLines))
        return(list(
            Total_Deduplicated_Percentage = NA_real_,
            Sequence_Duplication_Levels = data.frame(NULL)
        ))

    if (length(fqcLines[[mod]]) == 0)
        return(list(
            Total_Deduplicated_Percentage = NA_real_,
            Sequence_Duplication_Levels = data.frame(NULL)
        ))

    x <- fqcLines[[mod]]

    ## Check for the presence of the Total Deduplicate Percentage value
    hasTotDeDup <- grepl("Total Deduplicated Percentage", x)
    Total_Deduplicated_Percentage <- NA_real_
    if (any(hasTotDeDup)) {
        Total_Deduplicated_Percentage <-
            gsub(".+\\t(.*)", "\\1", x[hasTotDeDup])
        Total_Deduplicated_Percentage <-
            as.numeric(Total_Deduplicated_Percentage)
    }
    if (length(Total_Deduplicated_Percentage) > 1)
        stop("Too many elements matched Total_Deduplicated_Percentage")

    ## Remove the Total value entry from the original object
    df <- .splitByTab(x[!hasTotDeDup])
    names(df) <- gsub(" ", "_", names(df))

    ## Check for the required values
    reqVals <- c(
        "Duplication_Level",
        ## "Percentage_of_deduplicated", # Omitted with FastQC > 0.12
        "Percentage_of_total"
    )
    cols <- names(df)
    stopifnot(reqVals %in% cols)

    ## Convert percentages to numeric
    num_cols <- setdiff(cols, "Duplication_Level")
    df[num_cols] <- lapply(df[num_cols], as.numeric)

    ## Return a list with both values
    list(
        Total_Deduplicated_Percentage = Total_Deduplicated_Percentage,
        Sequence_Duplication_Levels = as_tibble(df)
    )

}

.getOverrepSeq <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Overrepresented_sequences"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])
    names(df) <- gsub(" ", "_", names(df))

    ## Check for the required values
    reqVals <- c("Sequence", "Count", "Percentage", "Possible_Source")
    stopifnot(reqVals %in% names(df))

    df$Count <- as.integer(df$Count)
    df$Percentage <- as.numeric(df$Percentage)
    as_tibble(df)

}

.getAdapterCont <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Adapter_Content"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])
    names(df) <- gsub(" ", "_", names(df))

    ## Check for the required values
    reqVals <- "Position"
    stopifnot(reqVals %in% names(df))
    stopifnot(ncol(df) > 1)

    numCols <- !names(df) %in% reqVals
    df[numCols] <- lapply(df[numCols], as.numeric)
    as_tibble(df)

}

.getKmerCont <- function(fqcLines){

    ## Return empty df if the module is missing
    mod <- "Kmer_Content"
    if (!mod %in% names(fqcLines)) return(data.frame(NULL))
    if (length(fqcLines[[mod]]) == 0) return(data.frame(NULL))

    df <- .splitByTab(fqcLines[[mod]])
    names(df) <- gsub(" ", "_", names(df))

    ## Check for the required values
    reqVals <-
        c("Sequence", "Count", "PValue", "Obs/Exp_Max", "Max_Obs/Exp_Position")
    stopifnot(reqVals %in% names(df))

    df$Count <- as.integer(df$Count)
    df$PValue <- as.numeric(df$PValue)
    df$`Obs/Exp_Max` <- as.numeric(df$`Obs/Exp_Max`)

    as_tibble(df)

}
