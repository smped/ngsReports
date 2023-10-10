test_that(".FastpFile fails on incorrect structure", {
  fl <- system.file("extdata", "ATTG_R1_fastqc.zip", package = "ngsReports")
  expect_error(suppressWarnings(.FastpFile(fl)))
  expect_error(suppressWarnings(.FastpFile(NULL), "!is.null(x) is not TRUE"))
  expect_error(suppressWarnings(.FastpFile(""), "file.exists(x) is not TRUE"))

})

test_that("fastp file loads correctly", {
  fl <- system.file("extdata", "fastp.json.gz", package = "ngsReports")
  expect_true(is(.FastpFile(fl), ".FastpFile"))
  fp <- FastpData(fl)
  expect_true(is(fp, "FastpData"))

  summary <- getModule(fp, "Summary")
  expect_equal(
    names(summary), c("Before_filtering", "After_filtering", "Filtering_result")
  )
  expect_true(all(vapply(summary, is, logical(1), "data.frame")))

  adapters <- getModule(fp, "Adapters")
  expect_true(is(adapters, "data.frame"))
  expect_true(is(adapters$read1_adapter_count, "list"))

  dup <- getModule(fp, "Duplication")
  expect_true(is(dup, "data.frame"))
  expect_true(is(dup$histogram, "list"))

  ins <- getModule(fp, "Insert")
  expect_true(is(ins, "data.frame"))
  expect_true(is(ins$histogram, "list"))

  mod_cols <- c("quality_curves", "content_curves", "kmer_count", "overrepresented_sequences")
  b4 <- getModule(fp, "Before")
  expect_equal(names(b4), paste0("read", seq_len(2)))
  expect_true(all(vapply(b4, is, logical(1), "data.frame")))
  expect_true(all(vapply(b4$read1[mod_cols], is, logical(1), "list")))

  after <- getModule(fp, "After")
  expect_equal(names(after), paste0("read", seq_len(2)))
  expect_true(all(vapply(after, is, logical(1), "data.frame")))
  expect_true(all(vapply(after$read1[mod_cols], is, logical(1), "list")))

  expect_true(fp@paired)
  expect_equal(fp@version, NA_character_)

})

closeAllConnections()
