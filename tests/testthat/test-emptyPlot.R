context("Check structure of emptyPlot()")

p <- ngsReports:::emptyPlot("test")

test_that("plot can be printed", {
    expect_error(print(p), NA)
})

test_that("plot is ggplot",{
    expect_true(is.ggplot(p))
})

test_that("theme_void() is applied",{
    expect_equal(p$theme, theme_void())
})

test_that("data is empty",{
    expect_equal(p$data, structure(list(), class = "waiver"))
})

test_that("correct aesthetics",{
    l <- list(x = "x", y = "y", label = "x")
    expect_equal(p$labels, l)
})
