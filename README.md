[![Build Status](https://travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=Bioc3.6)](https://travis-ci.org/UofABioinformaticsHub/ngsReports)


# ngsReports

An R Package for managing FastQC reports and other NGS related log files inside R.
This branch is compatible with Bioconductor >= 3.7 only. To install this package using Bioconductor <= 3.6 (R <= 3.4.4) please use the drop down menu above to change to the branch Bioc3.6, and follow the instructions there.

## Installation
To install required packages follows the instructions below.
Currently you need to install the fastqcTheoreticalGC package separately.

**NB: This branch is currently frozen for Bioconductor <=3.6**.
To install using Bioconductor 3.7, select the branch `master` using the drop-down menu above, then follow the installation instructions on the corresponding page.

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocGenerics", "BiocStyle", "BSgenome", "checkmate", "devtools", "ggdendro",  "plotly", "reshape2", "Rsamtools", "scales", "ShortRead", "tidyverse",  "viridis", "viridisLite", "zoo"))
devtools::install_github('mikelove/fastqcTheoreticalGC')
devtools::install_github('UofABioinformaticsHub/ngsReports', ref = "Bioc3.6", build_vignettes = TRUE)

library(ngsReports)
```

# Vignette

The vignette for usage is [here](https://uofabioinformaticshub.github.io/ngsReports/vignettes/ngsReportsIntroduction)


