[![Build Status](https://api.travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=dev)](https://api.travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=dev)

# ngsReports

An R Package for managing FastQC reports and other NGS related log files inside R.
This branch is compatible with Bioconductor >= 3.7 only. To install this package using Bioconductor <= 3.6 (R <= 3.4.4) please use the drop down menu above to change to the branch Bioc3.6, and follow the instructions there.

## Installation
To install required packages follows the instructions below.
Currently you need to install the fastqcTheoreticalGC package separately.

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocGenerics", "BiocStyle", "BSgenome", "checkmate", "devtools", "ggdendro",  "plotly", "reshape2", "scales", "ShortRead", "tidyverse",  "viridisLite", "zoo"))
devtools::install_github('mikelove/fastqcTheoreticalGC')
devtools::install_github('UofABioinformaticsHub/ngsReports', "dev", build_vignettes = FALSE)
library(ngsReports)
```

# Vignette

The vignette for usage is [here](https://uofabioinformaticshub.github.io/ngsReports/vignettes/ngsReportsIntroduction)


