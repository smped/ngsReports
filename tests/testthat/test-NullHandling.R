context("Test correct handling of NULL objects")

test_that("FastqcData fails on NULL",{
    expect_error(FastqcData(NULL))
})

test_that("FastqcDataList fails on NULL",{
    expect_error(FastqcDataList(NULL))
})

test_that(".FastqcFile failes on NULL", {
    expect_error(.FastqcFile(NULL))
})