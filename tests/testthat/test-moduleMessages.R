context("Test that messages appear when modules are missing from some files")

f <- c(system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports"),
       system.file("extdata/errorTestingFiles/moduleNamesNoData.zip", package = "ngsReports")
)
fdl <- getFastqcData(f)

test_that("All modules give correct messages",{
  expect_message(Adapter_Content(fdl))
  expect_message(Basic_Statistics(fdl))
  expect_message(Kmer_Content(fdl))
  expect_message(Overrepresented_sequences(fdl))
  expect_message(Per_base_N_content(fdl))
  expect_message(Per_base_sequence_content(fdl))
  expect_message(Per_base_sequence_quality(fdl))
  expect_message(Per_sequence_GC_content(fdl))
  expect_message(Per_sequence_quality_scores(fdl))
  expect_message(Per_tile_sequence_quality(fdl))
  expect_message(Sequence_Duplication_Levels(fdl))
  expect_message(Sequence_Length_Distribution(fdl))
})

