context("Check correct behaviour for setLabels()")

test_that("setLabels errors correctly with missing Filename column",{
    df <- data.frame()
    expect_error(setLabels(df))
})

test_that("setLabels handles missing labels correctly",{
    df <- data.frame(Filename = "file.fq", stringsAsFactors = FALSE)
    expect_equal(setLabels(df), structure("file", names = "file.fq"))
})

test_that("setLabels handles label mismatches",{
    df <- data.frame(Filename = "", stringsAsFactors = FALSE)
    labs <- "x"
    expect_error(setLabels(df, labels = x))
})

test_that("setLabels handles duplicated labels",{
    df <- data.frame(Filename = rep("file.fq", 2))
    expect_error(setLabels(df))
})

test_that("setLabels correctly checks for data.frame",{
    df <- c()
    expect_error(setLabels(df))
})

test_that("setLabels correctly removes all suffixes",{
    suf <- c("fq", "fastq", "fq.gz", "fastq.gz", "bam", "sam", "cram")
    df <- data.frame(Filename = paste(LETTERS[seq_along(suf)], suf, sep = "."),
                     stringsAsFactors = FALSE)
    trueOut <- structure(LETTERS[seq_along(suf)], names = df$Filename)
    expect_equal(setLabels(df), trueOut)
})