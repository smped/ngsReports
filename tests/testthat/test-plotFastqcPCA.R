context("Check correct behaviour for plotFastqcPCA")

# Get the files included with the package
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)

# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)

test_that("group data.frame errors with incorrect columns",{
    groups <- data.frame()
    expect_error(
        plotFastqcPCA(fdl, module = "Per_sequence_quality_scores", cluster = TRUE, groups = groups),
        "c\\(\"Filename\", \"Group\"\\) %in% colnames\\(groups\\) are not all TRUE"
    )
})

test_that("group data.frame errors with missing values",{
    groups <- tibble(
        Filename = fqName(fdl),
        Group = stringr::str_extract(Filename, "R[12]")
    )[1,]
    expect_error(
        plotFastqcPCA(fdl, module = "Per_sequence_quality_scores", cluster = TRUE, groups = groups),
        "all\\(data\\$Filename %in% groups\\$Filename\\) is not TRUE"
    )
})

test_that("Missing module errors", {
    expect_error(plotFastqcPCA(fdl), "!missing\\(module\\) is not TRUE")
})

test_that("hulls works",{
    expect_s3_class(
        plotFastqcPCA(fdl, "Per_sequence_quality_scores", clusterType = "hulls"),
        c("gg", "ggplot")
    )
})

test_that("color works",{
    expect_s3_class(
        plotFastqcPCA(fdl, "Per_sequence_quality_scores", clusterType = "color"),
        c("gg", "ggplot")
    )
})
