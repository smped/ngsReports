# Immediate Issues

## Chris


## Hien


## Steve

- **Add methods for zero results for all plots, based on plotNContent()**
- Validation Functions for `fastqcData` and `fastqcDataLists` objects
- Build unit tests for classes and methods
- **Test on outlier FastQC reports (1 sequence etc)**
- Finish Vignette
    - Add summary tables
- **plotOverrepresentedSummary**
    - Check plot for no overrepresented sequences: *Temporary version added...*    
    - Include in default report
- Change plotly `plotReadTotals()` to use rectangles (as for `plotDuplicationLevels()`)
- Fix NOTES from `plotDuplicationLevels()` and write correctly
- Fix handling of zero OverrepresentedSequences in writeHTMLReport template file

# Important But Not Pressing

- Test the output of `runFastQC` now it has been changed to a `FastqcFileList`
- Set character methods for `runFastQC` & redefine this function as a method for class `FastqFileList`

# Future Plans

- Include the capacity for selective overwriting in `runFastQC`
- Calculate Ranks on each module and summarise
- Add a function to merge kMers
- Add Fastq Illumina Filter plot/status

## Chris
- PCA plot for >50 samples
- fix travis error
