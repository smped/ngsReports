test_that("Fails on file which is not fasta",{
    fl <- system.file("extdata", "bowtie2PE.txt", package = "ngsReports")
    expect_error(estGcDistn(fl))
})

faFile <- system.file("extdata", "Athaliana.TAIR10.tRNA.fasta", package = "ngsReports")

if (requireNamespace("truncnorm", quietly = TRUE)){

    test_that("Character method works as expected",{
        # This implicitly tests the DNAStringSet method
        df <- estGcDistn(
            faFile, n = 100, rl = 50, fl = 60, fragSd = 10, bins = 51
        )
        expect_equal(colnames(df), c("GC_Content", "Freq"))
        expect_equal(nrow(df), 51)
    })

    test_that("Key parameters error correctly",{
        expect_error(estGcDistn(faFile, n = 100, rl = 60, fl = 50))
        expect_error(estGcDistn(faFile, n = 100, rl = 50, fl = 60, fragSd = -1))
    })

}
