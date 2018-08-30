context("Test all functions for getting module data")

# Define the example file
x <- system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports")

test_that("Example file exists",{
  expect_true(file.exists(x))
})

# Extract the fastqc data
# Change this to test the lower-level functions.
# Read the lines in instead & test for the correct output
fqcData <- getFastqcData(x)

test_that("Basic Statistics is correct",{

  bs <- Basic_Statistics(fqcData)
  expect_true(
    setequal(names(bs), c("Filename", "Total_Sequences",
                          "Sequences_flagged_as_poor_quality",
                          "Shortest_sequence", "Longest_sequence",
                          "%GC", "File_type", "Encoding")) &&
      setequal(vapply(bs, typeof, character(1)),
               c("character", rep("integer", 4),rep("character", 3)))
  )
})

test_that("Per base sequence qualities is correct",{
  sq <- Per_base_sequence_quality(fqcData)
  expect_true(
    setequal(names(sq), c("Filename", "Base", "Mean", "Median",
                          "Lower_Quartile", "Upper_Quartile",
                          "10th_Percentile", "90th_Percentile")) &&
      setequal(vapply(sq,typeof, character(1)),
               c(rep("character", 2), rep("double", 6)))
  )
})
