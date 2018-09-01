context("Test all functions for correct errors and nulls")

# Define the example file with known structure
x <- system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports")

test_that("Example file exists",{
  expect_true(file.exists(x))
})

# Extract the fastqc data as the lines using the code from `getFastqcData()`
fl <- file.path( gsub(".zip$", "", basename(x)), "fastqc_data.txt")
uz <- unz(x, fl)
fastqcLines <- readLines(uz)
close(uz)

# Remove any '#' symbols
fastqcLines <- gsub("#", "", fastqcLines)
# Remove any lines which specify '>>END_MODULE'
fastqcLines <- fastqcLines[!grepl(">>END_MODULE", fastqcLines)]
# Remove the first (FastQC version) line
fastqcLines <- fastqcLines[-1]
# Split into modules
fastqcLines <- split(fastqcLines, cumsum(grepl("^>>", fastqcLines)))
# Assign the module names
names(fastqcLines) <- vapply(fastqcLines, function(x){
  # Get the first line before the tab separator
  nm <- gsub(">>(.+)\\t.+", "\\1", x[[1]])
  # Replace the space with underscore
  gsub(" ", "_", nm)
}, character(1))
# Remove the first line
fastqcLines <- lapply(fastqcLines, function(x){x[-1]})
# Now all the tests can be run

test_that("Example file is correct", {
  expect_true(
    setequal(names(fastqcLines),
             c("Basic_Statistics", "Per_base_sequence_quality",
               "Per_tile_sequence_quality", "Per_sequence_quality_scores",
               "Per_base_sequence_content", "Per_sequence_GC_content",
               "Per_base_N_content", "Sequence_Length_Distribution",
               "Sequence_Duplication_Levels", "Overrepresented_sequences",
               "Adapter_Content", "Kmer_Content"))
  )
})

test_that("Check getBasicStatistics() errors and nulls",{
  fastqcLines[["Basic_Statistics"]] <- fastqcLines[["Basic_Statistics"]][-1]
  otherMods <- names(fastqcLines) != "Basic_Statistics"
  expect_equal(getBasicStatistics(fastqcLines[otherMods]), data.frame(NULL))
  expect_error(getBasicStatistics(fastqcLines))
})

test_that("Check getPerBaseSeqQuals() errors and nulls",{
  mod <- "Per_base_sequence_quality"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getPerBaseSeqQuals(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getPerBaseSeqQuals(empty), data.frame(NULL))
  expect_error(getPerBaseSeqQuals(fastqcLines))
})

test_that("Check getPerTileSeqQuals() errors and nulls",{
  mod <- "Per_tile_sequence_quality"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getPerTileSeqQuals(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getPerTileSeqQuals(empty), data.frame(NULL))
  expect_error(getPerTileSeqQuals(fastqcLines))
})

test_that("Check getPerSeqQualScores() errors and nulls",{

  mod <- "Per_sequence_quality_scores"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getPerSeqQualScores(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getPerSeqQualScores(empty), data.frame(NULL))
  expect_error(getPerSeqQualScores(fastqcLines))
  
})

test_that("Check getPerBaseSeqContent() errors and nulls",{

  mod <- "Per_base_sequence_content"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getPerBaseSeqContent(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getPerBaseSeqContent(empty), data.frame(NULL))
  expect_error(getPerBaseSeqContent(fastqcLines))
  
})

test_that("Check getPerSeqGcContent() errors and nulls",{

  mod <- "Per_sequence_GC_content"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getPerSeqGcContent(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getPerSeqGcContent(empty), data.frame(NULL))
  expect_error(getPerSeqGcContent(fastqcLines))
  
})

test_that("Check getSeqLengthDist() errors and nulls",{

  mod <- "Sequence_Length_Distribution"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getSeqLengthDist(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getSeqLengthDist(empty), data.frame(NULL))
  expect_error(getSeqLengthDist(fastqcLines))
  
})

test_that("Check getKmerContent() errors and nulls",{

  mod <- "Kmer_Content"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getKmerContent(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getKmerContent(empty), data.frame(NULL))
  expect_error(getKmerContent(fastqcLines))
  
})

test_that("Check getAdapterContent() errors and nulls",{

  mod <- "Adapter_Content"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getAdapterContent(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getAdapterContent(empty), data.frame(NULL))
  expect_error(getAdapterContent(fastqcLines))
  
})

test_that("Check getOverrepSeq() errors and nulls",{
  mod <- "Overrepresented_sequences"
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(getOverrepSeq(fastqcLines[otherMods]), data.frame(NULL))
  expect_equal(getOverrepSeq(empty), data.frame(NULL))
  expect_error(getOverrepSeq(fastqcLines))
})

test_that("Check getSeqDuplicationLevels() errors and nulls",{
  
  mod <- "Sequence_Duplication_Levels"
  # Remove the first line
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  otherMods <- names(fastqcLines) != mod
  empty <- list(character(0))
  names(empty) <- mod
  expect_equal(
    getSeqDuplicationLevels(fastqcLines[otherMods]),
    list(Total_Deduplicated_Percentage = NA_real_, 
         Sequence_Duplication_Levels = data.frame(NULL))
  )
  expect_equal(
    getSeqDuplicationLevels(empty),
    list(Total_Deduplicated_Percentage = NA_real_, 
         Sequence_Duplication_Levels = data.frame(NULL))
  )
  expect_true(
    is.na(getSeqDuplicationLevels(fastqcLines)[["Total_Deduplicated_Percentage"]])
  )
  # Now remove the second line
  fastqcLines[[mod]] <- fastqcLines[[mod]][-1]
  expect_error(getSeqDuplicationLevels(fastqcLines))
})

test_that("Completely empty fastqc file errors",{
  f <- system.file("extdata/errorTestingFiles/completelyEmpty.zip", package = "ngsReports")
  fqcFile <- FastqcFile(f)
  expect_error(getFastqcData(fqcFile))
  expect_error(getSummary(fqcFile))
})