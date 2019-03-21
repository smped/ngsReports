context("Check all import functions behave and error correctly")

bowtieLogs <- system.file("extdata", c("bowtiePE.txt", "bowtieSE.txt"), package = "ngsReports")
bowtie2Logs <- system.file("extdata", c("bowtie2PE.txt", "bowtie2SE.txt"), package = "ngsReports")
dupLogs <- system.file("extdata", "Sample1_Dedup_metrics.txt", package = "ngsReports")
starLog <- system.file("extdata", "log.final.out", package = "ngsReports")

test_that("importBowtieLogs loads correctly",{
    df <- importBowtieLogs(bowtieLogs)
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
    df <- importHisat2Logs(bowtie2Logs)
    nm <- c("Filename", "Total_Reads", "Paired_Reads", "Unique_Unpaired",
            "Unique_In_Pairs", "Unique_Discordant_Pairs", "Multiple_Unpaired",
            "Multiple_In_Pairs", "Not_Aligned", "Alignment_Rate")
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(bowtie2Logs), df$Filename)
})

test_that("importStarLogs loads correctly",{
    ## NB: This is identical to importBowtie2Logs
    df <- importStarLogs(starLog)
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

test_that("importDupMetrics loads correctly",{
    dup <- importDupMetrics(dupLogs)
    nm <- c("Library", "Unpaired Reads Examined", "Read Pairs Examined",
            "Secondary Or Supplementary Rds", "Unmapped Reads", "Unpaired Read Duplicates",
            "Read Pair Duplicates", "Read Pair Optical Duplicates", "Percent Duplication",
            "Estimated Library Size")
    expect_equal(names(dup), c("metrics", "histogram"))
    expect_equal(colnames(dup$histogram), c("Library", "Bin", "Value"))
    expect_equal(colnames(dup$metrics), nm)
})

test_that("importBowtieLogs errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importBowtieLogs(bowtie2Logs))
})

test_that("importHisat2Logs errors correctly",{
    ## These are bowtie logs so should error
    expect_error(importHisat2Logs(bowtieLogs))
})

test_that("importStarLogs errors correctly",{
    ## These are bowtie logs so should error
    expect_error(importStarLogs(bowtieLogs))
})

test_that("importDupMetrics errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importDupMetrics(bowtie2Logs))
})
