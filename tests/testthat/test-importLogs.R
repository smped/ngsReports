context("Check all import functions behave and error correctly")

bowtieLogs <- system.file("extdata", c("bowtiePE.txt", "bowtieSE.txt"), package = "ngsReports")
bowtie2Logs <- system.file("extdata", c("bowtie2PE.txt", "bowtie2SE.txt"), package = "ngsReports")
dupLogs <- system.file("extdata", "Sample1_Dedup_metrics.txt", package = "ngsReports")
starLog <- system.file("extdata", "log.final.out", package = "ngsReports")
arFile <- system.file("extdata", "adapterRemoval.settings", package = "ngsReports")
fcFile <- system.file("extdata", "featureCounts.summary", package = "ngsReports")
caFiles <- system.file("extdata", c("cutadapt_full.txt", "cutadapt_minimal.txt"), package = "ngsReports")

test_that("importNgsLogs fails on empty",{
    expect_error(importNgsLogs(bowtieLogs[0], "bowtie"))
})

test_that("importBowtieLogs works correctly",{
    df <- importNgsLogs(bowtieLogs, type = "bowtie")
    nm <- c("Filename", "Reads_Processed",
            "Reads_With_At_Least_One_Reported_Alignment",
            "Reads_That_Failed_To_Align",
            "Reads_With_Alignments_Suppressed_Due_To_-m",
            "Time_Loading_Reference",
            "Time_Loading_Forward_Index",
            "Time_Loading_Mirror_Index",
            "Time_Searching",
            "Overall_Time",
            "Seeded_Quality_Full_Index_Search"
            )
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(bowtieLogs), df$Filename)
})

test_that("importHisat2Logs loads correctly",{
    ## NB: This is identical to importBowtie2Logs
    df <- importNgsLogs(bowtie2Logs, type = "bowtie2")
    nm <- c("Filename", "Total_Reads", "Paired_Reads", "Unique_Unpaired",
            "Unique_In_Pairs", "Unique_Discordant_Pairs", "Multiple_Unpaired",
            "Multiple_In_Pairs", "Not_Aligned", "Alignment_Rate")
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(bowtie2Logs), df$Filename)
})

test_that("importStarLogs loads correctly",{
    ## NB: This is identical to importBowtie2Logs
    df <- importNgsLogs(starLog, type = "star")
    nm <- c("Filename", "Total_Mapped_Percent", "Number_Of_Input_Reads",
            "Average_Input_Read_Length", "Uniquely_Mapped_Reads_Number",
            "Uniquely_Mapped_Reads_Percent", "Average_Mapped_Length", "Number_Of_Reads_Mapped_To_Multiple_Loci",
            "Percent_Of_Reads_Mapped_To_Multiple_Loci", "Number_Of_Reads_Mapped_To_Too_Many_Loci",
            "Percent_Of_Reads_Mapped_To_Too_Many_Loci", "Percent_Of_Reads_Unmapped_Too_Many_Mismatches",
            "Percent_Of_Reads_Unmapped_Too_Short", "Percent_Of_Reads_Unmapped_Other",
            "Number_Of_Splices_Total", "Number_Of_Splices_Annotated_Sjdb",
            "Number_Of_Splices_GT/AG", "Number_Of_Splices_GC/AG", "Number_Of_Splices_AT/AC",
            "Number_Of_Splices_Non-Canonical", "Started_Job_On", "Started_Mapping_On",
            "Finished_On", "Mapping_Duration", "Mapping_Speed_Million_Of_Reads_Per_Hour",
            "Mismatch_Rate_Per_Base_Percent", "Deletion_Rate_Per_Base", "Deletion_Average_Length",
            "Insertion_Rate_Per_Base", "Insertion_Average_Length")
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(starLog), df$Filename)
})

nm <- c("LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED",
        "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS",
        "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE")
test_that("importDupMetrics loads correctly",{

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", 1)
    expect_equal(colnames(dup), nm)

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", 2)
    expect_equal(colnames(dup), c("LIBRARY", "BIN", "VALUE"))
})

test_that("importDupMetrics handles which as a character",{

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", "metrics")
    expect_equal(colnames(dup), nm)

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", "hist")
    expect_equal(colnames(dup), c("LIBRARY", "BIN", "VALUE"))
})

test_that("importBowtieLogs errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importNgsLogs(bowtie2Logs, type = "bowtie"))
})

test_that("importHisat2Logs errors correctly",{
    ## These are bowtie logs so should error
    expect_error(importNgsLogs(bowtieLogs, type = "hisat2"))
})

test_that("importStarLogs errors correctly",{
    ## These are bowtie logs so should error
    expect_error(importNgsLogs(bowtieLogs, type = "star"))
})

test_that("importDupMetrics errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importNgsLogs(bowtie2Logs, type = "duplicationMetrics"))
})

test_that("importAdapterRemovalLogs errors correctly",{
    expect_error(importNgsLogs(bowtie2Logs, type = "adapterRemoval"))
})

test_that("parseAdapterRemovalLogs parses correctly", {
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "settings")), 1L)
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "sequence")), 1L)
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "statistics")), 1L)
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "distribution")), 77L)
})

test_that("importFeatureCountsLogs errors correctly", {
    expect_error(importNgsLogs(bowtie2Logs, type = "featureCounts"))
})

test_that("parseAdapterRemovalLogs parses correctly", {
    cols <- c(
        "Filename", "Sample", "Assigned", "Unassigned_Ambiguity",
        "Unassigned_MultiMapping", "Unassigned_NoFeatures",
        "Unassigned_Unmapped", "Unassigned_MappingQuality",
        "Unassigned_FragmentLength", "Unassigned_Chimera",
        "Unassigned_Secondary", "Unassigned_Nonjunction",
        "Unassigned_Duplicate")
    expect_equal(nrow(importNgsLogs(fcFile, "feature")), 8L)
    expect_equal(ncol(importNgsLogs(fcFile, "feature")), 13L)
    expect_equal(colnames(importNgsLogs(fcFile, "feature")), cols)
})

test_that("importCutadapt errors correctly", {
    expect_error(importNgsLogs(bowtie2Logs, "cutadapt"))
    expect_error(importNgsLogs(caFiles, "cutadapt"))
    expect_error(importNgsLogs(caFiles[[1]], "cutadapt", which = 3))
    expect_error(importNgsLogs(caFiles[[1]], "cutadapt", which = 4))
    expect_message(importNgsLogs(caFiles[[2]], "cutadapt", which = 2))
})

test_that("parseCutadaptLogs parses correctly",{
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = 1)), c(1, 9)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = "summary")), c(1, 9)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = 2)), c(1, 8)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = "adapter1")), c(1, 8)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = 5)), c(49, 6)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = "overview")), c(49, 6)
    )

})
