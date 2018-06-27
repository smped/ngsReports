# Immediate Issues

- Tests on all hidden functions
- Get build time down
- Resolve all issues from `R CMD BiocCheck`

## Chris
- **Test on outlier FastQC reports (1 sequence etc)**
- Start commenting all functions better so Steve can understand them

## Hien


## Steve

- Build unit tests for classes and methods
- Finish Vignette
    - Add summary tables
- Fix `scale_fill_pwf()`


# Important But Not Pressing

- Test the output of `runFastQC` now it has been changed to a `FastqcFileList`
- Set character methods for `runFastQC` & redefine this function as a method for class `FastqFileList`

# Future Plans

- Include the capacity for selective overwriting in `runFastQC`
- Calculate Ranks on each module and summarise
- Add a function to merge kMers
- Add Fastq Illumina Filter plot/status
- BigData shiny app (PCA?)

