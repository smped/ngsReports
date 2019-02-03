context("Check Theoretical GC Content functions")

test_that("Fails on file which is not fasta",{
    fl <- system.file("extdata", "bowtie2PE.txt", package = "ngsReports")
    expect_error(getGcDistribution(fl))
})

faFile <- system.file("extdata", "Athaliana.TAIR10.tRNA.fasta", package = "ngsReports")

test_that("Character method works as expected",{
    df <- getGcDistribution(faFile, n = 100, readLength = 50, fragLength = 60, fragSd = 10, bins = 51)
    expect_equal(colnames(df), c("GC_Content", "Freq"))
    expect_equal(nrow(df), 51)
})

test_that("DNAStringSet method works as expected",{
    ss <- readDNAStringSet(faFile)
    df <- getGcDistribution(faFile, n = 100, readLength = 50, fragLength = 60, fragSd = 10, bins = 51)
    expect_equal(colnames(df), c("GC_Content", "Freq"))
    expect_equal(nrow(df), 51)
})

test_that("Key parameters error correctly",{
    expect_error(getGcDistribution(faFile, n = 100, readLength = 60, fragLength = 50))
    expect_error(getGcDistribution(faFile, n = 100, readLength = 50, fragLength = 60, fragSd = 0))
})
