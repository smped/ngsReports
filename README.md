
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
## Usage 
The paper for the package can be found [here](https://doi.org/10.1093/bioinformatics/btz937). 
For a detailed usage guide please see [here](https://bioconductor.org/packages/release/bioc/vignettes/ngsReports/inst/doc/ngsReportsIntroduction.html).


## ShinyApp

A Graphical User Interface (Shiny App) has been developed for interactive inspection of many FastQC reports. The ngsReports shiny app can be installed [here](https://github.com/UofABioinformaticsHub/shinyNgsReports).

## Tools Supported (By Category)

### Quality control
- FastQC
### Adapter removal and trimming
- AdapterRemoval
- CutAdapt
- Trimmomatic
### Mapping and alignment 
- Bowtie
- Bowtie2
- HISAT2
- Samtools flagstat
- STAR
- Picard MarkDuplicates
### Transcript/gene quantificaiton
- Feature counts
### Genome assembly
- BUSCO
- Quast

# Citation 

Please cite our [paper](https://doi.org/10.1093/bioinformatics/btz937):

```
@article{ward2019ngsreports,
  title={ngsReports: A Bioconductor package for managing FastQC reports and other NGS related log files},
  author={Ward, Christopher M and Thu-Hien, To and Pederson, Stephen M},
  journal={Bioinformatics}
}

```
