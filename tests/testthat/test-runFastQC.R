context("Check .fastqc()")

## Load a Fastq File
library(ShortRead)
sp <- SolexaPath(system.file('extdata', package = 'ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
f <- FastqFile(fl)

test_that("runFastqc Errors Correctly on missing executable",{
    expect_error(runFastQC(f, tempdir(), exec = ""))
    expect_error(runFastQC(f, tempdir(), exec = NULL))
})

test_that("runFastQC handles incorrect file types before searching for the executable",{
    logFile <- system.file("extdata", "log.final.out", package = "ngsReports")
    expect_error(out <- runFastQC(logFile, tempdir()))
})

## Check the function runs if the executable exists
exec <- Sys.which("fastqc")
if (exec != "") {
    f <- FastqFile(fl)
    ff <- suppressMessages(runFastQC(f, tempdir(), exec = exec))
    expect_true(methods::is(ff, "FastqcFile"))
}

## This ensures that any stray connections formed during testing are closed
closeAllConnections()