
[![Build Status](https://travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=master)](https://travis-ci.org/UofABioinformaticsHub/ngsReports.svg?branch=master)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)


# ngsReports

An R Package for managing FastQC reports and other NGS related log files inside R.
Except for some periodic minor bug fixes, this branch is the current release which is also available from [the Bioconductor website](https://bioconductor.org/packages/release/bioc/html/ngsReports.html).
To install this package for Bioconductor <= 3.6 (R <= 3.4.4) please use the drop down menu above to change to the branch Bioc3.6, and follow the instructions there.

## Installation

To install required packages follow the instructions below.

```
install.packages("BiocManager")
BiocManager::install("UofABioinformaticsHub/ngsReports")
library(ngsReports)
```

## ShinyApp

A Graphical User Interface (Shiny App) has been developed for interactive inspection of many FastQC reports. The ngsReports shiny app can be installed [here](https://github.com/UofABioinformaticsHub/shinyNgsReports).

# Citation 

Please cite our [preprint](https://www.biorxiv.org/content/early/2018/05/02/313148):

```
@article{ward2018ngsreports,
  title={ngsReports: An R Package for managing FastQC reports and other NGS related log files.},
  author={Ward, Christopher M and To, Hien and Pederson, Stephen M},
  journal={bioRxiv},
  pages={313148},
  year={2018},
  publisher={Cold Spring Harbor Laboratory}
}
```
