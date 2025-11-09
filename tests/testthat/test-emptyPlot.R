p <- ngsReports:::.emptyPlot("test")

test_that("plot can be printed", {
    expect_error(print(p), NA)
})

test_that("plot is ggplot",{
    expect_true(is_ggplot(p))
})

test_that("theme_void() is applied",{
    expect_equal(p@theme, theme_void())
})

test_that("data is empty",{
    expect_equal(p@data, structure(list(), class = "waiver"))
    expect_equal(
        p@layers$geom_text$aes_params$label, "test"
    )
})

test_that("Empty labels",{
    expect_true(length(p@labels) == 0)
})
