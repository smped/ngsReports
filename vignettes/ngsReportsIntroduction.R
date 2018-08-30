## --------------------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE,warning = FALSE)

## --------------------------------------------------------------------------
library(ngsReports)

## ---- eval = FALSE---------------------------------------------------------
#  fileDir <- system.file("extdata", package = "ngsReports")
#  writeHtmlReport(fileDir)

## ---- eval = FALSE---------------------------------------------------------
#  altTemplate <- file.path("path", "to", "template.Rmd")
#  writeHtmlReport(fileDir, template = altTemplate)

## --------------------------------------------------------------------------
fileDir <- system.file("extdata", package = "ngsReports")
files <- list.files(fileDir, pattern = "fastqc.zip$", full.names = TRUE)
fdl <- getFastqcData(files)

## ---- results='hide'-------------------------------------------------------
reads <- readTotals(fdl)

## --------------------------------------------------------------------------
library(dplyr)
library(pander)
filter(reads, grepl("R1", Filename)) %>% 
  pander(big.mark = ",")

## ----plotSummary, fig.cap="Default summary of FastQC flags.", fig.wide = TRUE----
plotSummary(fdl)

## --------------------------------------------------------------------------
plotReadTotals(fdl)

## --------------------------------------------------------------------------
plotReadTotals(fdl) +
  geom_vline(xintercept = 25000, linetype = 2) 

## --------------------------------------------------------------------------
plotBaseQualities(fdl)

## --------------------------------------------------------------------------
plotBaseQualities(fdl[[1]])

## --------------------------------------------------------------------------
plotBaseQualities(fdl[1:4], plotType = "boxplot")

## --------------------------------------------------------------------------
plotSequenceQualities(fdl)

## --------------------------------------------------------------------------
r2 <- seq(2, 12, by = 2)
plotSequenceQualities(fdl[r2], plotType = "line")

## --------------------------------------------------------------------------
plotSequenceContent(fdl)

## --------------------------------------------------------------------------
plotSequenceContent(fdl[[1]])

## --------------------------------------------------------------------------
plotNContent(fdl)

## --------------------------------------------------------------------------
plotSequenceLengthDistribution(fdl)

## --------------------------------------------------------------------------
plotDuplicationLevels(fdl)

## --------------------------------------------------------------------------
plotAdapterContent(fdl)

## --------------------------------------------------------------------------
plotAdapterContent(fdl[[1]])

## --------------------------------------------------------------------------
plotKmers(fdl)

## --------------------------------------------------------------------------
plotKmers(fdl[[1]], axis.text.x = element_text(angle = 90)) 

## ---- eval=FALSE-----------------------------------------------------------
#  transcriptomes(gcTheoretical)
#  genomes(gcTheoretical)

## --------------------------------------------------------------------------
plotGcContent(fdl)

## --------------------------------------------------------------------------
plotGcContent(fdl, theoreticalType = "Transcriptome", species = "Mmusculus")

## ----message=FALSE,warning=FALSE-------------------------------------------
faFile <- system.file("extdata", "Athaliana.TAIR10.tRNA.fasta", package="ngsReports")
plotGcContent(fdl, Fastafile = faFile)

## --------------------------------------------------------------------------
plotGcContent(fdl, plotType = "line",  theoreticalType = "Transcriptome")

## --------------------------------------------------------------------------
plotOverrepresentedSummary(fdl)

## ---- fig.wide = TRUE------------------------------------------------------
plotOverrepresentedSummary(fdl[[1]])

## ---- eval = FALSE---------------------------------------------------------
#  exportOverrepresented(fdl, n = 10)

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

