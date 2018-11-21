context("Check correct behaviour for modifying PwfCols objects")

pwf <- ngsReports::pwf
pwfWithAlpha <- setAlpha(pwf, 0.5)
vals <- c("#00CC00", "#E6E633", "#CC3333", "#FFFFFF")
names(vals) <- c("PASS", "WARN", "FAIL", "MAX")

test_that("getColours() works correctly",{
    expect_equal(getColours(pwf), vals)   
})

test_that("setAlpha() works correctly",{
    expect_equal(getColours(pwfWithAlpha), 
                 structure(paste0(vals, "80"), names = names(vals)))
})

test_that("Set Alpha ignores existing alpha values",{
    expect_equal(setAlpha(pwf, 0.6), setAlpha(pwfWithAlpha, 0.6))
})

test_that("setAlpha rejects invalid alpha values",{
    expect_error(setAlpha(pwf, -1))
    expect_error(setAlpha(pwf, "a"))
})

test_that("setColours rejects non-RGB colours",{
    expect_warning(setColours(pwf, PASS = "green"))
    expect_error(setColours(pwf, PASS = 1))
})

