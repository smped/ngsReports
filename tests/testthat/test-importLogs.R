context("Check all import functions behave and error correctly")

bowtieLogs <- system.file("extdata", c("bowtiePE.log", "bowtieSE.log"), package = "ngsReports")
bowtie2Logs <- system.file("extdata", c("bowtie2PE.log", "bowtie2SE.log"), package = "ngsReports")

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

test_that("importBowtieLogs errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importBowtieLogs(bowtie2Logs))
})
