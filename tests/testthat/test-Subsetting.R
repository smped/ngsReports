packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1:2]
fdl <- FastqcDataList(fl)

test_that("Single bracket subsetting works",{
    expect_length(fdl[1], 1)
    expect_length(fdl[names(fdl)[1]], 1)
})

test_that("Double bracket subsetting works",{
    expect_true(is(fdl[[1]], "FastqcData"))
    expect_true(is(fdl[[names(fdl)[1]]], "FastqcData"))
})

test_that("Subsetting out of range errors",{
    expect_error(fdl[0])
    expect_error(fdl[3])
    expect_error(fdl["foo"])
})
