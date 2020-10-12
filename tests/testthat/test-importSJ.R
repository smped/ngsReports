sjFiles <- system.file("extdata", "SJ.out.tab", package = "ngsReports")
df <- importSJ(sjFiles)

test_that(
    "Correct dimensions of import",
    expect_equal(dim(df), c(10, 11))
)

test_that(
    "Error on NULL", expect_error(importSJ(c()))
)

bwt <- system.file("extdata", "bowtiePE.txt", package = "ngsReports")
test_that(
    "Error on incorrect format",
    expect_error(importSJ(bwt))
)