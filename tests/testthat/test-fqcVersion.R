context("Test fqcVersion")

packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1:2]

# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)

test_that("fqcVersion gives correct output",{
    expect_equal(fqcVersion(fdl[[1]]), "0.11.2")
    expect_equal(fqcVersion(fdl)$version, rep("0.11.2", 2))
})

test_that("fqcVersion handles unimplements classes",{
    expect_message(
        fqcVersion(c()), "Method not implemented for objects of class NULL"
    )
})