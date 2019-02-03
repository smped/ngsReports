## ---- echo=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## --------------------------------------------------------------------------
library(ngsReports)

## --------------------------------------------------------------------------
fileDir <- system.file("extdata", package = "ngsReports")
files <- list.files(fileDir, pattern = "fastqc.zip$", full.names = TRUE)
fdl <- getFastqcData(files)

## ---- results='hide'-------------------------------------------------------
reads <- readTotals(fdl)

## --------------------------------------------------------------------------
library(dplyr)
library(pander)
dplyr::filter(reads, grepl("R1", Filename)) %>% 
  pander(big.mark = ",", caption = "Read totals from R1 libraries", justify = "lr")

## ----plotSummary, fig.cap="Default summary of FastQC flags.", fig.wide = TRUE----
plotSummary(fdl)

## --------------------------------------------------------------------------
plotReadTotals(fdl)

## --------------------------------------------------------------------------
plotReadTotals(fdl) +
    theme(
        legend.position = c(1, 0), 
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))

## --------------------------------------------------------------------------
plotBaseQualities(fdl)

## --------------------------------------------------------------------------
plotBaseQualities(fdl[[1]])

## --------------------------------------------------------------------------
plotBaseQualities(fdl[1:4], plotType = "boxplot")

## --------------------------------------------------------------------------
plotSequenceQualities(fdl)

## --------------------------------------------------------------------------
r2 <- grepl("R2", fileName(fdl))
plotSequenceQualities(fdl[r2], plotType = "line")

## --------------------------------------------------------------------------
plotSequenceContent(fdl)

## --------------------------------------------------------------------------
plotSequenceContent(fdl[[1]])

## --------------------------------------------------------------------------
plotSequenceContent(fdl[1:2], plotType = "line", nc = 1)

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

## ---- message=FALSE, warning=FALSE-----------------------------------------
faFile <- system.file("extdata", "Athaliana.TAIR10.tRNA.fasta", package = "ngsReports")
plotGcContent(fdl, Fastafile = faFile, n = 1000)

## --------------------------------------------------------------------------
plotGcContent(fdl, plotType = "line",  theoreticalType = "Transcriptome")

## --------------------------------------------------------------------------
plotOverrepresentedSummary(fdl)

## ---- fig.wide = TRUE------------------------------------------------------
plotOverrepresentedSummary(fdl[[1]])

## ---- eval = FALSE---------------------------------------------------------
#  exportOverrepresented(fdl, n = 10)

## --------------------------------------------------------------------------
fl <- c("bowtie2PE.txt", "bowtie2SE.txt")
bowtie2Logs <- system.file("extdata", fl, package = "ngsReports")
df <- importBowtie2Logs(bowtie2Logs)

## ---- echo=FALSE-----------------------------------------------------------
pander(df, style = "rmarkdown", caption = "Example of SE and PE output from bowtie 2")

## --------------------------------------------------------------------------
fls <- c("bowtiePE.txt", "bowtieSE.txt")
bowtieLogs <- system.file("extdata", fls, package = "ngsReports")
df <- importBowtieLogs(bowtieLogs)

## --------------------------------------------------------------------------
starLog <- system.file("extdata", "log.final.out", package = "ngsReports")
df <- importStarLogs(starLog)

## --------------------------------------------------------------------------
sysDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(sysDir, "Dedup_metrics.txt", full.names = TRUE)
dupMetrics <- importDuplicationMetrics(fl)
str(dupMetrics)

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

