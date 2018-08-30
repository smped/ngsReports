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
  expect_null(getBasicStatistics(fastqcLines[otherMods]))
  expect_error(getBasicStatistics(fastqcLines))
})

test_that("Check getPerBaseSeqQuals() errors and nulls",{
  fastqcLines[["Per_base_sequence_quality"]] <- fastqcLines[["Per_base_sequence_quality"]][-1]
  otherMods <- names(fastqcLines) != "Per_base_sequence_quality"
  expect_null(getPerBaseSeqQuals(fastqcLines[otherMods]))
  expect_error(getPerBaseSeqQuals(fastqcLines))
})

test_that("Check getPerTileSeqQuals() errors and nulls",{
  fastqcLines[["Per_tile_sequence_quality"]] <- fastqcLines[["Per_tile_sequence_quality"]][-1]
  otherMods <- names(fastqcLines) != "Per_tile_sequence_quality"
  expect_null(getPerTileSeqQuals(fastqcLines[otherMods]))
  expect_error(getPerTileSeqQuals(fastqcLines))
})

test_that("Check getPerSeqQualScores() errors and nulls",{
  fastqcLines[["Per_sequence_quality_scores"]] <- fastqcLines[["Per_sequence_quality_scores"]][-1]
  otherMods <- names(fastqcLines) != "Per_sequence_quality_scores"
  expect_null(getPerSeqQualScores(fastqcLines[otherMods]))
  expect_error(getPerSeqQualScores(fastqcLines))
})

test_that("Check getPerBaseSeqContent() errors and nulls",{
  fastqcLines[["Per_base_sequence_content"]] <- fastqcLines[["Per_base_sequence_content"]][-1]
  otherMods <- names(fastqcLines) != "Per_base_sequence_content"
  expect_null(getPerBaseSeqContent(fastqcLines[otherMods]))
  expect_error(getPerBaseSeqContent(fastqcLines))
})

test_that("Check getPerSeqGcContent() errors and nulls",{
  fastqcLines[["Per_sequence_GC_content"]] <- fastqcLines[["Per_sequence_GC_content"]][-1]
  otherMods <- names(fastqcLines) != "Per_sequence_GC_content"
  expect_null(getPerSeqGcContent(fastqcLines[otherMods]))
  expect_error(getPerSeqGcContent(fastqcLines))
})

test_that("Check getSeqLengthDist() errors and nulls",{
  fastqcLines[["Sequence_Length_Distribution"]] <- fastqcLines[["Sequence_Length_Distribution"]][-1]
  otherMods <- names(fastqcLines) != "Sequence_Length_Distribution"
  expect_null(getSeqLengthDist(fastqcLines[otherMods]))
  expect_error(getSeqLengthDist(fastqcLines))
})

test_that("Check getKmerContent() errors and nulls",{
  fastqcLines[["Kmer_Content"]] <- fastqcLines[["Kmer_Content"]][-1]
  otherMods <- names(fastqcLines) != "Kmer_Content"
  expect_null(getKmerContent(fastqcLines[otherMods]))
  expect_error(getKmerContent(fastqcLines))
})

test_that("Check getAdapterContent() errors and nulls",{
  fastqcLines[["Adapter_Content"]] <- fastqcLines[["Adapter_Content"]][-1]
  otherMods <- names(fastqcLines) != "Adapter_Content"
  expect_null(getAdapterContent(fastqcLines[otherMods]))
  expect_error(getAdapterContent(fastqcLines))
})