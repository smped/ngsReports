## This just perfoms generic tests for the correct objects for each function
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
# Load the FASTQC data as a FastqcDataList object
fdl <- FastqcDataList(fl)
fd <- FastqcData(fl[[1]])

fp <- system.file("extdata", "fastp.json.gz", package = "ngsReports") |>
  FastpData()
fpl <- FastpDataList(path(fp))
