x <- c(
    list.files("inst/extdata/", pattern = "samt", full.names = TRUE)
)

## Need to sort out the filename & bamfilename yet
.parseSamtoolsStatsLogs <- function(data, which = 1){

    allFiles <- lapply(data, function(x){

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
    })

    ## Check all  modules
    modNames <- lapply(allFiles, names)
    modNames <- unique(unlist(modNames))

    out <- lapply(modNames, function(x){
        df <- lapply(allFiles, function(f){
            if (x %in% names(f)) f[[x]]
            else tibble(Filename = unique(f[["CHK"]][["Filename"]]))
        })
        dplyr::bind_rows(df)
    })
    names(out) <- modNames
    out
}

#' @title Check for correct structure of supplied samtools stats files
#' @description Check for correct structure of supplied samtools stats
#' files after reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.isValidSamtoolsStatsLog <- function(x){
    chkLn1 <- grepl("This file was produced by samtools stats", x[[1]])
    chkSN <- any(grepl("^# Summary Numbers.", x))
    all(chkLn1, chkSN)
}

## Try figuring out how to get a boxplot of quality scores
getBxp <- function(df){
    df <- mutate(df, Q = cumsum(Count) / sum(Count))
    Q25 <- filter(df, Q > 0.25)$Score[1]
    Q50 <- filter(df, Q > 0.5)$Score[1]
    Q75 <- filter(df, Q > 0.75)$Score[1]
    iqr <- Q75 - Q25
    lwr <- max(min(df$Score), Q25 - 1.5*iqr)
    upr <- min(max(df$Score), Q75 + 1.5*iqr)
    mn <- with(df, sum(Score * Count) / sum(Count))
    df[1:2][1,] %>%
        mutate(lwr = lwr, Q25 = Q25, Q50 = Q50, Q75 = Q75, upr = upr, mn = mn)
}
