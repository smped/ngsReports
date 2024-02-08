test_that(".FastqcFile can be formed from zip archive", {
  fl <- system.file("extdata", "ATTG_R1_fastqc.zip", package = "ngsReports")
  fqcFile <- .FastqcFile(fl)
  expect_equal(length(fqcFile), 1)
})

test_that(".FastqcFile can be formed from directory",{
    fl <- system.file("extdata", "exampleRNASeq_R1_fastqc", package = "ngsReports")
    fqcFile <- .FastqcFile(fl)
    expect_equal(length(fqcFile), 1)
})

test_that(".FastqcFile fails on incorrect directory structure", {
    d <- system.file("extdata", package = "ngsReports")
    expect_error(suppressWarnings(.FastqcFile(d)))
})

closeAllConnections()
