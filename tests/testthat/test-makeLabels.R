context("Check correct behaviour for makeLabels()")

test_that("makeLabels errors correctly with missing Filename column",{
    df <- data.frame()
    expect_error(makeLabels(df))
})

test_that("makeLabels handles missing labels correctly",{
    df <- data.frame(Filename = "file.fq", stringsAsFactors = FALSE)
    expect_equal(makeLabels(df), structure("file", names = "file.fq"))
})

test_that("makeLabels handles label mismatches",{
    df <- data.frame(Filename = "", stringsAsFactors = FALSE)
    labs <- "x"
    expect_error(makeLabels(df, labels = labs))
})

test_that("makeLabels handles duplicated labels",{
    df <- data.frame(Filename = rep("file.fq", 2))
    expect_error(makeLabels(df))
})

test_that("makeLabels correctly checks for data.frame",{
    df <- c()
    expect_error(makeLabels(df))
})

test_that("makeLabels correctly removes all suffixes",{
    suf <- c("fq", "fastq", "fq.gz", "fastq.gz", "bam", "sam", "cram")
    df <- data.frame(Filename = paste(LETTERS[seq_along(suf)], suf, sep = "."),
                     stringsAsFactors = FALSE)
    trueOut <- structure(LETTERS[seq_along(suf)], names = df$Filename)
    expect_equal(makeLabels(df), trueOut)
})