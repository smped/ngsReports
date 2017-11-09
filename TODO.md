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
 
- Once the S4 methods are up (see below), add the individual `plotKmer` plots
- Add individual plots for `Sequence Duplication Levels`
- Fix `Error: attempt to select less than one element in get1index` in `Per Base Sequence Content` on first click on the page
    - Change base colours to be consistent between individual plots & the heatmap
- Fix `Per Sequence Quality Scores`: **Error**: Column `Filename` is unknown on first click   
- Change slider label to 'Max Sequences' for `Overrepresented Sequences`

## Individual Functions

- **Adapter Content**
    - add dendrogram when `usePlotly == TRUE`
- **plotKmers**
    - Fix plotly `NA` vals
    - Migrate to S4
    - Fix legend
- **plotNContent**
    - Migrate to S4
    - Remove `subset` argument
    - Merge `plotNContentPlotly()` with main function
    - Check plot for completely missing N content
- **plotOverrepresentedHeatmap**
    - Migrate to S4
    - Remove `subset` argument
    - Merge `plotOverrepresentedHeatmapPlotly()` with main function
    - Check plot for no overrepresented sequences
- **plotOverrepresentedSummary**
    - Migrate to S4
    - Remove `subset` argument
    - Merge `plotOverrepresentedHeatmapPlotly()` with main function
    - Check plot for no overrepresented sequences    
    - Include in default report & shiny app
- **plotSequenceContent**
    - Migrate to S4
    - Add clustering & dendrograms to plot AND shiny app
- **plotSummary**
    - Add clustering/dendrogram & tidy up for the app
