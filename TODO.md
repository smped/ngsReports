- Include the capacity for selective overwriting in `runFastQC`
- Set character methods for `runFastQC` & redefine this function as a method for class `FastqFileList`
- Test the output of `runFastQC` now it hs been changed to a `FastqcFileList`
- Validation Functions for `fastqcData` and `fastqcDataLists` objects
- Calculate Ranks on each module and summarise
- Add a function to merge kMers
- Add Fastq Illumina Filter plot/status
- Have a better look at http://multiqc.info/examples/rna-seq/multiqc_report.html#
- Build tests for classes and methods
- Finish Vignette

## Shiny App
 
- Fix `Error: attempt to select less than one element in get1index` in `Per Base Sequence Content` on first click on the page
    - Change base colours to be consistent between individual plots & the heatmap
- Fix `Per Sequence Quality Scores`: **Error**: Column `Filename` is unknown on first click   
- Set Overrepresented Summary to produce a table on click *added plot for now*
- Add `exportOverrepresented()`

## Individual Functions

- **plotOverrepresentedSummary**
    - Check plot for no overrepresented sequences: *Temporary version added...*    
    - Include in default report
- **plotSequenceContent**
    - Add clustering & dendrograms to plot AND shiny app
    
**Add methods for zero results for all plots, based on plotNContent()**
