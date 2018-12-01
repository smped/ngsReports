context("Test all functions for correct module data")

# Define the example file with known structure
x <- system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports")

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

test_that("Check getBasicStatistics() is correct",{
  bs <- getBasicStatistics(fqcLines)
  expect_true(
    setequal(names(bs),
             c("Filename", "Total_Sequences", "Sequences_flagged_as_poor_quality",
               "Shortest_sequence", "Longest_sequence", "%GC", "File_type",
               "Encoding")))
  expect_true(
    all(vapply(bs, typeof, character(1)) == c(
      "character", rep("integer", 4), rep("character", 3)
    )))
  expect_equal(nrow(bs), 1)
})

test_that("Check getPerBaseSeqQuals() is correct",{
  sq <- getPerBaseSeqQuals(fqcLines)
  expect_true(
    setequal(
      names(sq),
      c("Base", "Mean", "Median", "Lower_Quartile", "Upper_Quartile",
        "10th_Percentile","90th_Percentile")
      ))
  expect_true(
      all(
        vapply(sq, typeof, character(1)) == c("character", rep("double", 6))
        ))
  expect_equal(nrow(sq), 47)
})

test_that("Check getPerTileSeqQuals() is correct",{
  sq <- getPerTileSeqQuals(fqcLines)
  expect_true(
    setequal(
      names(sq), c("Tile", "Base", "Mean")
    ))
  expect_true(
    all(
      vapply(sq, typeof, character(1)) == c("character", "character", "double")
    ))
  expect_equal(nrow(sq), 4512)
})

test_that("Check getPerSeqQualScores() is correct",{
  sq <- getPerSeqQualScores(fqcLines)
  expect_true(
    setequal(names(sq), c("Quality","Count"))
    )
  expect_true(
    all(vapply(sq, is.integer, logical(1)))
    )
  expect_equal(nrow(sq), 39)
})

test_that("Check getPerBaseSeqCont() is correct",{
  sc <- getPerBaseSeqCont(fqcLines)
  expect_true(
    setequal(names(sc), c("Base", "G", "A", "T", "C"))
  )
  expect_true(
    all(vapply(sc, typeof, character(1)) == c("character", rep("double", 4)))
  )
  expect_equal(nrow(sc), 47)
})

test_that("Check getPerSeqGcCont() is correct",{
  gc <- getPerSeqGcCont(fqcLines)
  expect_true(
    setequal(names(gc), c("GC_Content", "Count"))
  )
  expect_true(
    all(vapply(gc, typeof, character(1)) == c("integer", "double"))
  )
  expect_equal(nrow(gc), 101)
})

test_that("Check getSeqLengthDist() is correct",{
  df <- getSeqLengthDist(fqcLines)
  expect_true(
    setequal(names(df), c("Length", "Lower", "Upper", "Count"))
  )
  expect_true(
    all(vapply(df, typeof, character(1)) == c("character", rep("integer", 3)))
  )
  expect_equal(nrow(df), 1)
})

test_that("Check getSeqDuplicationLevels() provides correct output",{
  res <- getSeqDuplicationLevels(fqcLines)
  expect_true(
    setequal(names(res),
             c("Total_Deduplicated_Percentage", "Sequence_Duplication_Levels"))
  )
  expect_equal(typeof(res[["Total_Deduplicated_Percentage"]]), "double")
  expect_true(
    setequal(
      names(res[["Sequence_Duplication_Levels"]]),
      c("Duplication_Level", "Percentage_of_deduplicated", "Percentage_of_total")
    )
  )
  expect_true(
    all(
      vapply(res[["Sequence_Duplication_Levels"]], typeof, character(1)) == c(
        "character", "double", "double"
      )
    )
  )
  expect_true(
    setequal(res[["Sequence_Duplication_Levels"]][["Duplication_Level"]],
             c(1:9, paste0(">", c(10, 50, 100, 500, "1k", "5k", "10k+"))))
  )

})

test_that("Check getOverrepSeq() is correct",{
  df <- getOverrepSeq(fqcLines)
  expect_true(
    setequal(names(df), c("Sequence", "Count", "Percentage", "Possible_Source"))
  )
  expect_true(
    all(
      vapply(df, typeof, character(1)) == c("character", "integer", "double", "character")
    )
  )
  expect_equal(nrow(df), 182)
})


test_that("Check getAdapterCont() is correct",{
  df <- getAdapterCont(fqcLines)
  nAdapTypes <- ncol(df) - 1
  expect_equal(names(df)[1], "Position")
  expect_gt(ncol(df), 1)
  expect_true(
    all(
      vapply(df, typeof, character(1)) == c("character",
                                            rep("double", nAdapTypes))
    )
  )
  expect_equal(nrow(df), 72)
})

test_that("Check getKmerCont() is correct",{
  df <- getKmerCont(fqcLines)
  expect_true(
    setequal(names(df),
             c("Sequence", "Count", "PValue", "Obs/Exp_Max",
               "Max_Obs/Exp_Position")
    )
  )
  expect_true(
    all(
      vapply(df, typeof, character(1)) == c("character", "integer",
                                            "double", "double", "character")
    )
  )
  expect_equal(nrow(df), 20)
})
