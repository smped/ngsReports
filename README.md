[![Build Status](https://api.travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=dev)](https://api.travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=dev)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)


# ngsReports

An R Package for managing FastQC reports and other NGS related log files inside R.
This branch is compatible with Bioconductor >= 3.7 only. To install this package using Bioconductor <= 3.6 (R <= 3.4.4) please use the drop down menu above to change to the branch Bioc3.6, and follow the instructions there.

## Installation
To install required packages follows the instructions below.

```
install.packages("BiocManager")
pkgs <- c(
    "BiocGenerics", "BiocStyle", "BSgenome", "checkmate", "ggdendro", "plotly", 
    "reshape2", "scales", "ShortRead", "tidyverse", "viridisLite", "zoo", 
    "mikelove/fastqcTheoreticalGC"
)
BiocManager::install(pkgs)
BiocManager::install("UofABioinformaticsHub/ngsReports", ref = "dev")
library(ngsReports)
```

# Vignette

The vignette for usage is [here](https://uofabioinformaticshub.github.io/ngsReports/vignettes/ngsReportsIntroduction)


