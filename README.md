
![Build Status](https://github.com/steveped/ngsReports/workflows/check-bioc/badge.svg)
![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)


# ngsReports

An R Package for managing FastQC reports and other NGS related log files inside R.
Except for some periodic minor bug fixes, this branch is the current release which is also available from [the Bioconductor website](https://bioconductor.org/packages/release/bioc/html/ngsReports.html).
To install this package for Bioconductor <= 3.6 (R <= 3.4.4) please use the drop down menu above to change to the branch Bioc3.6, and follow the instructions there.
Versions for Bioc releases 3.9 and 3.10 are also available as branches from this repository.

## Installation

To install required packages follow the instructions below.

```
install.packages("BiocManager")
BiocManager::install("steveped/ngsReports")
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
@article{ward2018ngsreports,
    author = {Ward, Christopher M and To, Thu-Hien and Pederson, Stephen M},
    title = "{ngsReports: a Bioconductor package for managing FastQC reports and other NGS related log files}",
    journal = {Bioinformatics},
    year = {2019},
    month = {12},
    doi = {10.1093/bioinformatics/btz937},
    url = {https://doi.org/10.1093/bioinformatics/btz937}
}
```
