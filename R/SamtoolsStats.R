#' @title The SamtoolsStats Object Class
#'
#' @description The SamtoolsStats Object Class
#'
#' @details This object class is the main object required for generating plots
#' and tables. Modules are
#' contained as individual slots, which can be viewed using \code{slotNames}.
#'
#' Individual modules can be returned using the function \code{getModule()}
#' and specifying which module is required. See \code{\link{getModule}} for
#' more details.
#'
#' All modules contain tibbles, except for cmd
#'
#' @slot Summary_Numbers Summary of alignments
#' @slot First_Fragment_Qualities Distribution of quality scores by cycle
#' @slot Last_Fragment_Qualities Distribution of quality scores by cycle
#' @slot GC_Content_of_first_fragments As above but only for first fragments
#' @slot GC_Content_of_last_fragments As above but only for last fragments
#' @slot ACGT_content_per_cycle Base content per cycle
#' @slot ACGT_content_for_first_fragments As above but only for first fragments
#' @slot ACGT_content_for_last_fragments As above but only for last fragments
#' @slot Insert_sizes Breakdown of insert sizes and alignments
#' @slot Read_lengths Read length distribution
#' @slot Read_lengths_first_fragments As above but only for first fragments
#' @slot Read_lengths_last_fragments As above but only for last fragments
#' @slot Indels_per_cycle Forward and Reverse insertions/deletions per cycle
#' @slot Coverage_distribution Coverage by bin
#' @slot GC_depth GC distribution with 10th, 25th, 50th, 75th and 90th
#' percentiles
#' @slot Checksums CRC32 Checksums for Read Names, Sequences & Qualities
#' @slot cmd The command used for generation of the report
#' @slot version The version of samtools & htslib used for generation of the
#' report
#'
#' @return An object of class SamtoolsStats
#'
#' @include validationFunctions.R
#'
#' @rdname SamtoolsStats
#' @aliases SamtoolsStats-class
setClass(
    "SamtoolsStats",
    slots = c(
        Summary_Numbers = "data.frame",
        First_Fragment_Qualities = "data.frame",
        Last_Fragment_Qualities = "data.frame",
        GC_Content_of_first_fragments = "data.frame",
        GC_Content_of_last_fragments = "data.frame",
        ACGT_content_per_cycle = "data.frame",
        ACGT_content_for_first_fragments = "data.frame",
        ACGT_content_for_last_fragments = "data.frame",
        Insert_sizes = "data.frame",
        Read_lengths = "data.frame",
        Read_lengths_first_fragments = "data.frame",
        Read_lengths_last_fragments = "data.frame",
        Indels_per_cycle = "data.frame",
        Coverage_distribution = "data.frame",
        GC_depth = "data.frame",
        Checksums = "data.frame",
        cmd = "character",
        version = "data.frame"
    )
)
setValidity("SamtoolsStats", .isValidSamtoolsStats)

#' @param x Path to a single samtools stats report
#' @rdname SamtoolsStats
#' @export
SamtoolsStats <- function(x){
    ## Some basic Checks
    stopifnot(!is.null(x))
    stopifnot(length(x) == 1)
    stopifnot(file.exists(x))

    ## Now parse the data
    ssLines <- readLines(x)
    stopifnot(.isValidSamtoolsStatsFile(x))
    ssList <- .parseSamtoolsStats(x)
}

.isValidSamtoolsStatsFile <- function(x){
    ## Check for the presence of all modules in parsed lines
    checks <- c(
        L1 = grepl("This file was produced by samtools stats", x[[1]]),
        SN = any(grepl("^# Summary Numbers.", x)),
        FFQ = any(grepl("^# First Fragment Qualities.", x)),
        LFQ = any(grepl("^# Last Fragment Qualities", x)),
        GCF = any(grepl("^# GC Content of first fragments", x)),
        GCL = any(grepl("^# GC Content of last fragments", x)),
        GCC = any(grepl("^# ACGT content per cycle", x)),
        FBC = any(grepl("^# ACGT content .+ first frag", x)),
        LBC = any(grepl("^# ACGT content .+ last frag", x)),
        IS = any(grepl("^# Insert sizes", x)),
        RL = any(grepl("^# Read lengths.", x)),
        FRL = any(grepl("^# Read lengths - first", x)),
        LRL = any(grepl("^# Read lengths - last", x)),
        ID = any(grepl("^# Indel distribution", x)),
        IC = any(grepl("^# Indels per cycle", x)),
        COV = any(grepl("^# Coverage distribution", x)),
        GCD = any(grepl("^# GC-depth", x))
    )
    all(checks)
}

## Need to sort out the filename & bamfilename yet
.parseSamtoolsStats <- function(x){

    ## Get the bam/sam Filename
    bamName <- stringr::str_subset(x, "The command line was")
    bamName <- stringr::str_remove(bamName, ".+The command line was: +stats ")

    ## Define the modules contained inthe file
    modLns <- stringr::str_subset(x, "Use `grep")
    mods <- stringr::str_replace(modLns, ".+ Use `grep \\^([^ ]+) .+", "\\1")
    mods <- c("CHK", mods)

    ## Split into modules, then remove any lines from each which are not
    ## tab delimited
    splitLines <- split(x, f = cumsum(stringr::str_detect(x, "Use `grep")))
    names(splitLines) <- mods
    splitLines[-1] <- lapply(splitLines[-1], function(x){
        ## Just keep the lines with a tab delimiter
        x <- stringr::str_subset(x, "\t")
        ## Remove any trailing comments
        x <- stringr::str_remove(x, "\t#.+")
        ## Remove the first column
        .splitByTab(x, FALSE)[-1]
    })

    ## The Summary Numbers module
    if ("SN" %in% mods) {
        nm <- stringr::str_remove(splitLines[["SN"]]$V2, ":$")
        nm <- stringr::str_replace_all(nm, " ", "_")
        splitLines[["SN"]] <- structure(
            as.numeric(splitLines[["SN"]]$V3),
            names = nm
        )
        splitLines[["SN"]] <- as_tibble(as.list(splitLines[["SN"]]))
        splitLines[["SN"]] <- dplyr::mutate_at(
            splitLines[["SN"]], vars(starts_with("is_")), as.logical
        )
        splitLines[["SN"]] <- dplyr::mutate_at(
            splitLines[["SN"]],
            vars(
                starts_with("reads_"),
                ends_with("sequences"),
                ends_with("fragments"),
                starts_with("pairs"),
                ends_with("pairs"),
                contains("non-primary")
            ),
            as.integer
        )
        splitLines[["SN"]]$Filename <- bamName
        splitLines[["SN"]] <- dplyr::select(
            splitLines[["SN"]], "Filename", tidyselect::everything()
        )
    }

    ## First Fragment Qualities
    if ("FFQ" %in% mods) {
        colnames(splitLines[["FFQ"]]) <- c(
            "cycle", paste0("Q", seq_len(ncol(splitLines[["FFQ"]]) - 1))
        )
        splitLines[["FFQ"]] <- as_tibble(splitLines[["FFQ"]])
        splitLines[["FFQ"]] <-
            dplyr::mutate_all(splitLines[["FFQ"]], as.integer)
        splitLines[["FFQ"]]$Filename <- bamName
        splitLines[["FFQ"]] <- dplyr::select(
            splitLines[["FFQ"]], "Filename", tidyselect::everything()
        )
        splitLines[["FFQ"]] <- tidyr::gather(
            splitLines[["FFQ"]], "score", "count", starts_with("Q")
        )
        splitLines[["FFQ"]][["score"]] <- str_remove(
            splitLines[["FFQ"]][["score"]], "Q"
        )
        splitLines[["FFQ"]][["score"]] <-
            as.integer(splitLines[["FFQ"]][["score"]])
        splitLines[["FFQ"]] <-
            dplyr::filter(splitLines[["FFQ"]], count > 0)
    }

    ## Last Fragment Qualities
    if ("LFQ" %in% mods) {
        colnames(splitLines[["LFQ"]]) <- c(
            "cycle", paste0("Q", seq_len(ncol(splitLines[["LFQ"]]) - 1))
        )
        splitLines[["LFQ"]] <- as_tibble(splitLines[["LFQ"]])
        splitLines[["LFQ"]] <-
            dplyr::mutate_all(splitLines[["LFQ"]], as.integer)
        splitLines[["LFQ"]]$Filename <- bamName
        splitLines[["LFQ"]] <- dplyr::select(
            splitLines[["LFQ"]], "Filename", tidyselect::everything()
        )
        splitLines[["LFQ"]] <- tidyr::gather(
            splitLines[["LFQ"]], "score", "count", starts_with("Q")
        )
        splitLines[["LFQ"]][["score"]] <- str_remove(
            splitLines[["LFQ"]][["score"]], "Q"
        )
        splitLines[["LFQ"]][["score"]] <-
            as.integer(splitLines[["LFQ"]][["score"]])
        splitLines[["LFQ"]] <-
            dplyr::filter(splitLines[["LFQ"]], count > 0)
    }

    ## GC Content of first fragments
    if ("GCF" %in% mods) {
        colnames(splitLines[["GCF"]]) <- c("gc", "count")
        splitLines[["GCF"]] <- as_tibble(splitLines[["GCF"]])
        splitLines[["GCF"]][["gc"]] <-
            as.numeric(splitLines[["GCF"]][["gc"]])
        splitLines[["GCF"]][["count"]] <-
            as.integer(splitLines[["GCF"]][["count"]])
        splitLines[["GCF"]]$Filename <- bamName
        splitLines[["GCF"]] <- dplyr::select(
            splitLines[["GCF"]], "Filename", tidyselect::everything()
        )

    }

    ## GC Content of last fragments
    if ("GCL" %in% mods) {
        colnames(splitLines[["GCL"]]) <- c("gc", "count")
        splitLines[["GCL"]] <- as_tibble(splitLines[["GCL"]])
        splitLines[["GCL"]][["gc"]] <-
            as.numeric(splitLines[["GCL"]][["gc"]])
        splitLines[["GCL"]][["count"]] <-
            as.integer(splitLines[["GCL"]][["count"]])
        splitLines[["GCL"]]$Filename <- bamName
        splitLines[["GCL"]] <- dplyr::select(
            splitLines[["GCL"]], "Filename", tidyselect::everything()
        )

    }

    ## ACGT content per cycle
    if ("GCC" %in% mods) {
        colnames(splitLines[["GCC"]]) <- c(
            "cycle", "A", "C", "G", "T", "N", "O"
        )
        splitLines[["GCC"]] <- as_tibble(splitLines[["GCC"]])
        splitLines[["GCC"]][["cycle"]] <-
            as.integer(splitLines[["GCC"]][["cycle"]])
        splitLines[["GCC"]] <- dplyr::mutate_at(
            splitLines[["GCC"]],
            c("A", "C", "G", "T", "N", "O"),
            as.numeric
        )
        splitLines[["GCC"]]$Filename <- bamName
        splitLines[["GCC"]] <- dplyr::select(
            splitLines[["GCC"]], "Filename", tidyselect::everything()
        )

    }

    ## ACGT content per cycle for first fragments
    if ("FBC" %in% mods) {
        colnames(splitLines[["FBC"]]) <- c(
            "cycle", "A", "C", "G", "T", "N", "O"
        )
        splitLines[["FBC"]] <- as_tibble(splitLines[["FBC"]])
        splitLines[["FBC"]][["cycle"]] <-
            as.integer(splitLines[["FBC"]][["cycle"]])
        splitLines[["FBC"]] <- dplyr::mutate_at(
            splitLines[["FBC"]],
            c("A", "C", "G", "T", "N", "O"),
            as.numeric
        )
        splitLines[["FBC"]]$Filename <- bamName
        splitLines[["FBC"]] <- dplyr::select(
            splitLines[["FBC"]], "Filename", tidyselect::everything()
        )

    }

    ## ACGT content per cycle for last fragments
    if ("LBC" %in% mods) {
        colnames(splitLines[["LBC"]]) <- c(
            "cycle", "A", "C", "G", "T", "N", "O"
        )
        splitLines[["LBC"]] <- as_tibble(splitLines[["LBC"]])
        splitLines[["LBC"]][["cycle"]] <-
            as.integer(splitLines[["LBC"]][["cycle"]])
        splitLines[["LBC"]] <- dplyr::mutate_at(
            splitLines[["LBC"]],
            c("A", "C", "G", "T", "N", "O"),
            as.numeric
        )
        splitLines[["LBC"]]$Filename <- bamName
        splitLines[["LBC"]] <- dplyr::select(
            splitLines[["LBC"]], "Filename", tidyselect::everything()
        )

    }

    ## Insert sizes
    if ("IS" %in% mods) {
        colnames(splitLines[["IS"]]) <- c(
            "insert_size", "pairs_total", "inward_oriented_pairs",
            "outward_oriented_pairs", "other_pairs"
        )
        splitLines[["IS"]] <- as_tibble(splitLines[["IS"]])
        splitLines[["IS"]] <- dplyr::mutate_all(
            splitLines[["IS"]], as.integer
        )
        splitLines[["IS"]]$Filename <- bamName
        splitLines[["IS"]] <- dplyr::select(
            splitLines[["IS"]], "Filename", tidyselect::everything()
        )
    }

    ## Read lengths
    if ("RL" %in% mods) {
        colnames(splitLines[["RL"]]) <- c("read_length", "count")
        splitLines[["RL"]] <- as_tibble(splitLines[["RL"]])
        splitLines[["RL"]] <- dplyr::mutate_all(
            splitLines[["RL"]], as.integer
        )
        splitLines[["RL"]]$Filename <- bamName
        splitLines[["RL"]] <- dplyr::select(
            splitLines[["RL"]], "Filename", tidyselect::everything()
        )

    }

    ## Read lengths - first fragments
    if ("FRL" %in% mods) {
        colnames(splitLines[["FRL"]]) <- c("read_length", "count")
        splitLines[["FRL"]] <- as_tibble(splitLines[["FRL"]])
        splitLines[["FRL"]] <- dplyr::mutate_all(
            splitLines[["FRL"]], as.integer
        )
        splitLines[["FRL"]]$Filename <- bamName
        splitLines[["FRL"]] <- dplyr::select(
            splitLines[["FRL"]], "Filename", tidyselect::everything()
        )

    }

    ## Read lengths - last fragments
    if ("LRL" %in% mods) {
        colnames(splitLines[["LRL"]]) <- c("read_length", "count")
        splitLines[["LRL"]] <- as_tibble(splitLines[["LRL"]])
        splitLines[["LRL"]] <- dplyr::mutate_all(
            splitLines[["LRL"]], as.integer
        )
        splitLines[["LRL"]]$Filename <- bamName
        splitLines[["LRL"]] <- dplyr::select(
            splitLines[["LRL"]], "Filename", tidyselect::everything()
        )

    }

    ## Indel distribution
    if ("ID" %in% mods) {
        colnames(splitLines[["ID"]]) <- c(
            "length", "number_of_insertions", "number_of_deletions"
        )
        splitLines[["ID"]] <- as_tibble(splitLines[["ID"]])
        splitLines[["ID"]] <- dplyr::mutate_all(
            splitLines[["ID"]], as.integer
        )
        splitLines[["ID"]]$Filename <- bamName
        splitLines[["ID"]] <- dplyr::select(
            splitLines[["ID"]], "Filename", tidyselect::everything()
        )
    }

    ## Indels per cycle
    if ("IC" %in% mods) {
        colnames(splitLines[["IC"]]) <- c(
            "cycle", "number_of_insertions_fwd", "number_of_insertions_rev",
            "number_of_deletions_fwd", "number_of_deletions_rev"
        )
        splitLines[["IC"]] <- as_tibble(splitLines[["IC"]])
        splitLines[["IC"]] <- dplyr::mutate_all(
            splitLines[["IC"]], as.integer
        )
        splitLines[["IC"]]$Filename <- bamName
        splitLines[["IC"]] <- dplyr::select(
            splitLines[["IC"]], "Filename", tidyselect::everything()
        )
    }

    ## Coverage distribution
    if ("COV" %in% mods) {
        colnames(splitLines[["COV"]]) <- c("range", "bin_max", "count")
        splitLines[["COV"]] <- as_tibble(splitLines[["COV"]])
        splitLines[["COV"]] <- dplyr::mutate_at(
            splitLines[["COV"]], c("bin_max", "count"), as.integer
        )
        splitLines[["COV"]]$Filename <- bamName
        splitLines[["COV"]] <- dplyr::select(
            splitLines[["COV"]], "Filename", tidyselect::everything()
        )
    }

    ## GC-depth
    if ("GCD" %in% mods) {
        colnames(splitLines[["GCD"]]) <- c(
            "gc", "unique_sequence_percentiles", "Q10", "Q25", "Q50", "Q75",
            "Q90"
        )
        splitLines[["GCD"]] <- as_tibble(splitLines[["GCD"]])
        splitLines[["GCD"]] <-
            dplyr::mutate_all(splitLines[["GCD"]], as.numeric)
        splitLines[["GCD"]]$Filename <- bamName
        splitLines[["GCD"]] <- dplyr::select(
            splitLines[["GCD"]], "Filename", tidyselect::everything()
        )
    }

    ## Now sort out the header section
    splitLines[["version"]] <- tibble(
        Filename = bamName,
        version = stringr::str_replace(
            splitLines[["CHK"]][[1]],
            ".*This file was produced by samtools stats \\(([^\\)]+)\\).+",
            "\\1"
        )
    )
    splitLines[["CHK"]] <- stringr::str_subset(
        splitLines[["CHK"]], "\t"
    )
    splitLines[["CHK"]] <- .splitByTab(splitLines[["CHK"]])
    splitLines[["CHK"]] <- splitLines[["CHK"]][-1]
    splitLines[["CHK"]] <- as_tibble(splitLines[["CHK"]])
    colnames(splitLines[["CHK"]]) <- stringr::str_remove(
        colnames(splitLines[["CHK"]]), "\\[[0-9]\\]"
    )
    colnames(splitLines[["CHK"]]) <- stringr::str_to_lower(
        colnames(splitLines[["CHK"]])
    )
    colnames(splitLines[["CHK"]]) <- stringr::str_replace_all(
        colnames(splitLines[["CHK"]]), " ", "_"
    )
    splitLines[["CHK"]]$Filename <- bamName
    splitLines[["CHK"]] <- dplyr::select(
        splitLines[["CHK"]], "Filename", tidyselect::everything()
    )

    splitLines

}
