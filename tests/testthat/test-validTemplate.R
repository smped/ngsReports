context("Check that supplied template is OK")

packageDir <- system.file("extdata", package = "ngsReports")
fileList <-
    list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1:6]

test_that("Supplied template compiles",{
    file.copy(fileList, tempdir())
    writeHtmlReport(tempdir(), usePlotly = FALSE, overwrite = TRUE)
    expect_true(file.exists(file.path(tempdir(), "ngsReports_Fastqc.html")))
})
