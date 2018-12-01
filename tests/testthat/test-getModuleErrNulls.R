context("Test all functions for correct errors and nulls")

# Define the example file with known structure
x <- system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports")

test_that("Example file exists",{
  expect_true(file.exists(x))
})

# Extract the fastqc data as the lines using the code from `getFastqcData()`
fl <- file.path( gsub(".zip$", "", basename(x)), "fastqc_data.txt")
uz <- unz(x, fl)
fqcLines <- readLines(uz)
close(uz)

# Remove any '#' symbols
fqcLines <- gsub("#", "", fqcLines)
# Remove any lines which specify '>>END_MODULE'
fqcLines <- fqcLines[!grepl(">>END_MODULE", fqcLines)]
# Remove the first (FastQC version) line
fqcLines <- fqcLines[-1]
# Split into modules
fqcLines <- split(fqcLines, cumsum(grepl("^>>", fqcLines)))
# Assign the module names
names(fqcLines) <- vapply(fqcLines, function(x){
  # Get the first line before the tab separator
  nm <- gsub(">>(.+)\\t.+", "\\1", x[[1]])
  # Replace the space with underscore
  gsub(" ", "_", nm)
}, character(1))
# Remove the first line
fqcLines <- lapply(fqcLines, function(x){x[-1]})
# Now all the tests can be run

test_that("Example file is correct", {
  expect_true(
    setequal(names(fqcLines),
             c("Basic_Statistics", "Per_base_sequence_quality",
               "Per_tile_sequence_quality", "Per_sequence_quality_scores",
               "Per_base_sequence_content", "Per_sequence_GC_content",
               "Per_base_N_content", "Sequence_Length_Distribution",
               "Sequence_Duplication_Levels", "Overrepresented_sequences",
               "Adapter_Content", "Kmer_Content"))
  )
})

test_that("Check .getBasicStatistics() errors and nulls",{
  fqcLines[["Basic_Statistics"]] <- fqcLines[["Basic_Statistics"]][-1]
  otherMods <- names(fqcLines) != "Basic_Statistics"
  expect_equal(.getBasicStatistics(fqcLines[otherMods]), data.frame(NULL))
  expect_error(.getBasicStatistics(fqcLines))
})

test_that("Check .getPerBaseSeqQuals() errors and nulls",{
  mod <- "Per_base_sequence_quality"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getPerBaseSeqQuals(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getPerBaseSeqQuals(empty), data.frame(NULL))
  expect_error(.getPerBaseSeqQuals(fqcLines))
})

test_that("Check .getPerTileSeqQuals() errors and nulls",{
  mod <- "Per_tile_sequence_quality"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getPerTileSeqQuals(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getPerTileSeqQuals(empty), data.frame(NULL))
  expect_error(.getPerTileSeqQuals(fqcLines))
})

test_that("Check .getPerSeqQualScores() errors and nulls",{

  mod <- "Per_sequence_quality_scores"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getPerSeqQualScores(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getPerSeqQualScores(empty), data.frame(NULL))
  expect_error(.getPerSeqQualScores(fqcLines))

})

test_that("Check .getPerBaseSeqCont() errors and nulls",{

  mod <- "Per_base_sequence_content"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getPerBaseSeqCont(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getPerBaseSeqCont(empty), data.frame(NULL))
  expect_error(.getPerBaseSeqCont(fqcLines))

})

test_that("Check .getPerSeqGcCont() errors and nulls",{

  mod <- "Per_sequence_GC_content"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getPerSeqGcCont(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getPerSeqGcCont(empty), data.frame(NULL))
  expect_error(.getPerSeqGcCont(fqcLines))

})

test_that("Check .getSeqLengthDist() errors and nulls",{

  mod <- "Sequence_Length_Distribution"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getSeqLengthDist(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getSeqLengthDist(empty), data.frame(NULL))
  expect_error(.getSeqLengthDist(fqcLines))

})

test_that("Check .getKmerCont() errors and nulls",{

  mod <- "Kmer_Content"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getKmerCont(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getKmerCont(empty), data.frame(NULL))
  expect_error(.getKmerCont(fqcLines))

})

test_that("Check .getAdapterCont() errors and nulls",{

  mod <- "Adapter_Content"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getAdapterCont(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getAdapterCont(empty), data.frame(NULL))
  expect_error(.getAdapterCont(fqcLines))

})

test_that("Check .getOverrepSeq() errors and nulls",{
  mod <- "Overrepresented_sequences"
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(.getOverrepSeq(fqcLines[otherMods]), data.frame(NULL))
  expect_equal(.getOverrepSeq(empty), data.frame(NULL))
  expect_error(.getOverrepSeq(fqcLines))
})

test_that("Check .getSeqDuplicationLevels() errors and nulls",{

  mod <- "Sequence_Duplication_Levels"
  # Remove the first line
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  otherMods <- names(fqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(
    .getSeqDuplicationLevels(fqcLines[otherMods]),
    list(Total_Deduplicated_Percentage = NA_real_,
         Sequence_Duplication_Levels = data.frame(NULL))
  )
  expect_equal(
    .getSeqDuplicationLevels(empty),
    list(Total_Deduplicated_Percentage = NA_real_,
         Sequence_Duplication_Levels = data.frame(NULL))
  )
  expect_true(
    is.na(.getSeqDuplicationLevels(fqcLines)[["Total_Deduplicated_Percentage"]])
  )
  # Now remove the second line
  fqcLines[[mod]] <- fqcLines[[mod]][-1]
  expect_error(.getSeqDuplicationLevels(fqcLines))
})

test_that("Completely empty fastqc file errors",{
  f <- system.file("extdata/errorTestingFiles/completelyEmpty.zip", package = "ngsReports")
  fqcFile <- FastqcFile(f)
  expect_error(getFastqcData(fqcFile))
  expect_error(.getSummary(fqcFile))
})

closeAllConnections()