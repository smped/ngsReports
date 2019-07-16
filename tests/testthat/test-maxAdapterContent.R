context("Check correct behaviour for maxAdapterContent()")

test_that("Error on object that is not a FastqcDataList",{
    expect_error(maxAdapterContent("a"))
})

