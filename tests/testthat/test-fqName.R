
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

test_that("Correct values returned for fastp files", {
    nm <- c(
        read1 = "S010_20170320003-4_ffpedna_pan-cancer-v1_S10_R1_001.fastq",
        read2 = "S010_20170320003-4_ffpedna_pan-cancer-v1_S10_R2_001.fastq"
    )
    expect_equal(fqName(fp), nm)
})
