fl <- system.file("extdata", "ATTG_R1_fastqc.zip", package = "ngsReports")
fd <- FastqcDataList(fl)

test_that("Correct values are returned",{
    expect_equal(fqName(fd), "ATTG_R1.fastq")
})

test_that("Incorrect values cannot be assigned",{
    expect_error(fqName(fd) <- 1)
    expect_error(fqName(fd) <- letters)
})

test_that("Correct values cannot be assigned",{
    fqName(fd) <- "a"
   expect_equal(fqName(fd), "a")
})
