[![Build Status](https://travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=master)](https://travis-ci.org/UofABioinformaticsHub/ngsReports)

# ngsReports

An R Package for managing FastQC reports and other NGS related log files inside R

## Installation
To install required packages 
`
packageList <- c(  "BiocGenerics", "checkmate", "dplyr", "ggdendro", "ggplot2",
"lubridate", "magrittr", "method"s, "plotly", "readr", "reshape2", "Rsamtools",
"scales", "shiny", "ShortRead", "stats", "stringr", "tibble", "viridis", "viridisLite",
"zoo", "shinyFiles")
packagesToInstall <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(packagesToInstall)) install.packages(packagesToInstall)
`

Install ngsReports from GitHub and load library
`
library(devtools)
install_github('UofABioinformaticsHub/ngsReports')
library(ngsReports)
`

# ShinyApp Usage 
For a analysis of multiple fastqc reports use the shinyApp by running:
`fastqcShiny()`
once inside the shiny app, files can be input by clicking the `Choose Files` button.
This will then open a pop-up window to select your fastqcReports (select multiple files by holding control, ect.) 
once selected files will load and first plot will appear.

ShinyApp can also be passed a character vector of filenames or a fastqcFileList
`
# Pass character vector fileList to shinyApp
 fileList <- list.files(path = "mydir/", pattern = ".zip$", full.names = TRUE)
 fastqcShiny(fastqcInput = fileList)
 
# Pass fastqcFileList to shinyApp
 fdl <- getFastqcData(fileList)
 fastqcShiny()
 
 `




