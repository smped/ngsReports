
test_that("Supplied template compiles",{
    file.copy(fl, tempdir(), overwrite = TRUE)
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
