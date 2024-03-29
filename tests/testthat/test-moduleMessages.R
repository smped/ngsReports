test_that("All modules give correct messages",{

  f <- c(
    system.file("extdata/ATTG_R1_fastqc.zip", package = "ngsReports"),
    system.file(
      "extdata/errorTestingFiles/moduleNamesNoData.zip",
      package = "ngsReports"
    )
  )
  fdl <- FastqcDataList(f)

  expect_message(getModule(fdl, "Adapter_Content"))
  expect_message(getModule(fdl, "Basic_Statistics"))
  expect_message(getModule(fdl, "Kmer_Content"))
  expect_message(getModule(fdl, "Overrepresented_sequences"))
  expect_message(getModule(fdl, "Per_base_N_content"))
  expect_message(getModule(fdl, "Per_base_sequence_content"))
  expect_message(getModule(fdl, "Per_base_sequence_quality"))
  expect_message(getModule(fdl, "Per_sequence_GC_content"))
  expect_message(getModule(fdl, "Per_sequence_quality_scores"))
  expect_message(getModule(fdl, "Per_tile_sequence_quality"))
  expect_message(getModule(fdl, "Sequence_Duplication_Levels"))
  expect_message(getModule(fdl, "Sequence_Length_Distribution"))
  expect_message(getModule(fdl, "Total_Deduplicated_Percentage"))
})

closeAllConnections()
