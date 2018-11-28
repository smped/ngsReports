# Immediate Issues

- Resolve all issues from `R CMD BiocCheck`

## Chris
- **Test on outlier FastQC reports (1 sequence etc)**
- Start commenting all functions better so Steve can understand them
- fix passing fdl to shiny app 
- Start developing Shiny app as it's own package

## Steve

- Check additional import functions
- Fix `scale_fill_pwf()`
- `maxAdapterContent()` as an S4 method


# Important But Not Pressing

- Test the output of `runFastQC` now it has been changed to a `FastqcFileList`
- Set character methods for `runFastQC` & redefine this function as a method for class `FastqFileList`

# Future Plans

- Include the capacity for selective overwriting in `runFastQC`
- Calculate Ranks on each module and summarise
- Add a function to merge kMers
- Add Fastq Illumina Filter plot/status
- BigData shiny app (PCA?)
- Add plots for imported log files
- Make padding using `plotly_empty()` optional for plots of a `FastqcData` object

