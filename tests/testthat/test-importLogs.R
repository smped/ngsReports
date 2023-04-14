bowtieLogs <- system.file("extdata", c("bowtiePE.txt", "bowtieSE.txt"), package = "ngsReports")
bowtie2Logs <- system.file("extdata", c("bowtie2PE.txt", "bowtie2SE.txt"), package = "ngsReports")
buscoFiles <- system.file("extdata", c("short_summary_Dmelanogaster_Busco.txt", "short_summary_Zcucurbitae_Busco.txt"), package = "ngsReports")
dupLogs <- system.file("extdata", "Sample1_Dedup_metrics.txt", package = "ngsReports")
starLog <- system.file("extdata", "log.final.out", package = "ngsReports")
quastFiles <- system.file("extdata", c("quast1.tsv", "quast2.tsv"), package = "ngsReports")
arFile <- system.file("extdata", "adapterRemoval.settings", package = "ngsReports")
fcFile <- system.file("extdata", "featureCounts.summary", package = "ngsReports")
caFiles <- system.file("extdata", c("cutadapt_full.txt", "cutadapt_minimal.txt"), package = "ngsReports")
trimoFiles <- system.file("extdata", "Sample1.trimmomaticPE.txt", package = "ngsReports")
flagFile <- system.file("extdata", "flagstat.txt", package = "ngsReports")
macs2File <- system.file("extdata", "macs2_callpeak.txt", package = "ngsReports")
umiFile <- system.file("extdata", "umitoolsDedup.txt", package = "ngsReports")

test_that("importNgsLogs fails on empty",{
    expect_error(importNgsLogs(bowtieLogs[0], "bowtie"))
})

test_that("importBowtieLogs works correctly",{
    df <- importNgsLogs(bowtieLogs, type = "bowtie")
    nm <- c("Filename", "Reads_Processed",
            "Reads_With_At_Least_One_Reported_Alignment",
            "Reads_That_Failed_To_Align",
            "Reads_With_Alignments_Suppressed_Due_To_-m",
            "Time_Loading_Reference",
            "Time_Loading_Forward_Index",
            "Time_Loading_Mirror_Index",
            "Time_Searching",
            "Overall_Time",
            "Seeded_Quality_Full_Index_Search"
            )
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(bowtieLogs), df$Filename)
})

test_that("importHisat2Logs loads correctly",{
    ## NB: This is identical to importBowtie2Logs
    df <- importNgsLogs(bowtie2Logs, type = "bowtie2")
    nm <- c("Filename", "Total_Reads", "Paired_Reads", "Unique_Unpaired",
            "Unique_In_Pairs", "Unique_Discordant_Pairs", "Multiple_Unpaired",
            "Multiple_In_Pairs", "Not_Aligned", "Alignment_Rate")
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(bowtie2Logs), df$Filename)
})

test_that("importStarLogs loads correctly",{
    ## NB: This is identical to importBowtie2Logs
    df <- importNgsLogs(starLog, type = "star")
    nm <- c("Filename", "Total_Mapped_Percent", "Number_Of_Input_Reads",
            "Average_Input_Read_Length", "Uniquely_Mapped_Reads_Number",
            "Uniquely_Mapped_Reads_Percent", "Average_Mapped_Length", "Number_Of_Reads_Mapped_To_Multiple_Loci",
            "Percent_Of_Reads_Mapped_To_Multiple_Loci", "Number_Of_Reads_Mapped_To_Too_Many_Loci",
            "Percent_Of_Reads_Mapped_To_Too_Many_Loci", "Percent_Of_Reads_Unmapped_Too_Many_Mismatches",
            "Percent_Of_Reads_Unmapped_Too_Short", "Percent_Of_Reads_Unmapped_Other",
            "Number_Of_Splices_Total", "Number_Of_Splices_Annotated_Sjdb",
            "Number_Of_Splices_GT/AG", "Number_Of_Splices_GC/AG", "Number_Of_Splices_AT/AC",
            "Number_Of_Splices_Non-Canonical", "Started_Job_On", "Started_Mapping_On",
            "Finished_On", "Mapping_Duration", "Mapping_Speed_Million_Of_Reads_Per_Hour",
            "Mismatch_Rate_Per_Base_Percent", "Deletion_Rate_Per_Base", "Deletion_Average_Length",
            "Insertion_Rate_Per_Base", "Insertion_Average_Length")
    ## Just check for the expected colnames & filenames
    expect_equal(colnames(df), nm)
    expect_equal(basename(starLog), df$Filename)
})


test_that("importQuastLog loads correctly", {
    df <- importNgsLogs(quastFiles, type = "quast")
    nm <- c("fileNames", "totalLength", "longestTig", "N50", "N75", "L50", "L75")
    expect_equal(colnames(df), nm)
})

test_that("importBuscoLog loads correctly", {
    df <- importNgsLogs(buscoFiles, type = "busco")
    nm <- c("name", "completeSingleCopy", "completeDuplicated", "fragmented", "missing")
    expect_equal(colnames(df), nm)
})

test_that("autodetect errors on multiple log types", {
    expect_error(importNgsLogs(c(bowtie2Logs, bowtieLogs), "auto"))
})

test_that("autodetect works", {
    x <- c(bowtieLogs, bowtie2Logs, buscoFiles, dupLogs, starLog, quastFiles, arFile, fcFile, caFiles)
    possTypes <- c(
        "adapterRemoval",
        "bowtie",
        "bowtie2",
        "busco",
        "cutadapt",
        "duplicationMetrics",
        "featureCounts",
        "hisat2",
        "quast",
        "star"
    )
    data <- suppressWarnings(lapply(x, readLines))
    names(data) <- basename(x)
    type <- vapply(data, .getToolName, character(1), possTypes = possTypes)
    tools <- c(
        rep("bowtie", length(bowtieLogs)),
        rep("bowtie2", length(bowtie2Logs)),
        rep("busco", length(buscoFiles)),
        rep("duplicationMetrics", length(dupLogs)),
        rep("star", length(starLog)),
        rep("quast", length(quastFiles)),
        rep("adapterRemoval", length(arFile)),
        rep("featureCounts", length(fcFile)),
        rep("cutadapt", length(caFiles))
    )

    expect_equal(as.character(type), tools)

    f <- system.file("extdata", "Athaliana.TAIR10.tRNA.fasta", package = "ngsReports")
    expect_error(importNgsLogs(f), "No matching file type was found")

})

nm <- c("LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED",
        "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS",
        "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE")
nm <- stringr::str_replace_all(nm, "_", " ")
nm <- stringr::str_to_title(nm)
test_that("importDupMetrics loads correctly",{

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", 1)
    expect_equal(colnames(dup), nm)

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", 2)
    expect_equal(colnames(dup), c("Library", "Bin", "Value"))
})

test_that("importDupMetrics handles which as a character",{

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", "metrics")
    expect_equal(colnames(dup), nm)

    dup <- importNgsLogs(dupLogs, "duplicationMetrics", "hist")
    expect_equal(colnames(dup), c("Library", "Bin", "Value"))
})

test_that("importBowtieLogs errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importNgsLogs(bowtie2Logs, type = "bowtie"))
})

test_that("importHisat2Logs errors correctly",{
    ## These are bowtie logs so should error
    expect_error(importNgsLogs(bowtieLogs, type = "hisat2"))
})

test_that("importStarLogs errors correctly",{
    ## These are bowtie logs so should error
    expect_error(importNgsLogs(bowtieLogs, type = "star"))
})

test_that("importDupMetrics errors correctly",{
    ## These are bowtie2 logs so should error
    expect_error(importNgsLogs(bowtie2Logs, type = "duplicationMetrics"))
})

test_that("importAdapterRemovalLogs errors correctly",{
    expect_error(importNgsLogs(bowtie2Logs, type = "adapterRemoval"))
})

test_that("parseAdapterRemovalLogs parses correctly", {
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "settings")), 1L)
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "sequence")), 1L)
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "statistics")), 1L)
    expect_equal(nrow(importNgsLogs(arFile, "adapter", "distribution")), 77L)
})

test_that("importFeatureCountsLogs errors correctly", {
    expect_error(importNgsLogs(bowtie2Logs, type = "featureCounts"))
})

test_that("parseAdapterRemovalLogs parses correctly", {
    cols <- c(
        "Filename", "Sample", "Assigned", "Unassigned_Ambiguity",
        "Unassigned_MultiMapping", "Unassigned_NoFeatures",
        "Unassigned_Unmapped", "Unassigned_MappingQuality",
        "Unassigned_FragmentLength", "Unassigned_Chimera",
        "Unassigned_Secondary", "Unassigned_Nonjunction",
        "Unassigned_Duplicate")
    expect_equal(nrow(importNgsLogs(fcFile, "feature")), 8L)
    expect_equal(ncol(importNgsLogs(fcFile, "feature")), 13L)
    expect_equal(colnames(importNgsLogs(fcFile, "feature")), cols)
})

test_that("importCutadapt errors correctly", {
    expect_error(importNgsLogs(bowtie2Logs, "cutadapt"))
    expect_error(importNgsLogs(caFiles, "cutadapt"))
    expect_error(importNgsLogs(caFiles[[1]], "cutadapt", which = 3))
    expect_message(importNgsLogs(caFiles[[2]], "cutadapt", which = 2))
})

test_that("parseCutadaptLogs parses correctly",{
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = 1)), c(1, 10)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = "summary")), c(1, 10)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = 2)), c(1, 10)
    )
    expect_equal(
        dim(importNgsLogs(caFiles[[1]], "cutadapt", which = "adapter1")), c(1, 10)
    )
})

test_that("parseTrimmomaticLogs behaves correctly",{
    df <- importNgsLogs(trimoFiles, "trimmomatic")
    expCols <- c(
        "Filename", "Type", "Input_Read_Pairs", "Both_Surviving",
        "Forward_Only_Surviving",  "Reverse_Only_Surviving", "Dropped",
        "Illumina_Clip", "Sliding_Window", "Trailing", "Min_Len",
        "Quality_Encoding"
    )
    expect_error(importNgsLogs(caFiles, "trimmomatic"))
    expect_equal(colnames(df), expCols)
})

test_that("parseFlagstatLogs behaves correctly", {

    df <- importNgsLogs(flagFile, "flagstat")
    expCols <- c("Filename", "QC-passed", "QC-failed", "flag")
    expect_equal(dim(df), c(13, 4))
    expect_equal(colnames(df), expCols)
})

test_that("parseMacs2CallPeak behaves correctly", {

    df <- importNgsLogs(macs2File, "macs2Callpeak")
    expCols <- c(
      "Filename", "name", "date", "paired_peaks", "min_length", "n_tags_treatment",
      "n_tags_control", "n_reads", "tag_length", "fragment_length",
      "alt_fragment_length", "format", "ChIP_seq_file", "control_file",
      "effective_genome_size", "band_width", "model_fold", "qvalue_cutoff",
      "max_gap", "keep_dup", "nomodel", "scale_to", "local", "broad",
      "paired_end", "outputs"
    )
    expect_equal(colnames(df), expCols)
    ## Parsing errors will appear as NA values
    expect_equal(
        sum(vapply(df, is.na, logical(1))), 0
    )

})

test_that("umitools dedup parses", {
  df <- importNgsLogs(umiFile, "umitoolsDedup")
  expCols <- c(
    "Filename", "Input_Reads", "Read_pairs", "Read_2_unmapped",
    "Number_of_reads_out", "Total_number_of_positions_deduplicated",
    "Mean_number_of_unique_UMIs_per_position", "Max_number_of_unique_UMIs_per_position",
    "Mates_never_found", "start", "duration"
  )
  expect_equal(colnames(df), expCols)
})
