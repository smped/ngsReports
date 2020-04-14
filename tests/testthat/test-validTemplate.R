context("Check that supplied template is OK")

packageDir <- system.file("extdata", package = "ngsReports")
fileList <-
    list.files(packageDir, pattern = "[ACGT]{4}_R[12].fastqc.zip", full.names = TRUE)[1:6]

test_that("Supplied template compiles",{
    file.copy(fileList, tempdir(), overwrite = TRUE)
    writeHtmlReport(tempdir(), usePlotly = FALSE, overwrite = TRUE)
    expect_true(file.exists(file.path(tempdir(), "ngsReports_Fastqc.html")))
})

test_that("writeHtmlReport Errors Correctly", {
    expect_error(writeHtmlReport(fastqcDir = ""))
    expect_error(writeHtmlReport(packageDir, template = ""))
    expect_error(
        writeHtmlReport(packageDir, file.path(packageDir, "hisat2PE.txt"))
    )
})
