# Immediate Issues

## Chris

- `Per_base_sequence_content` is selecting files in the reverse order
- `Overrepresented_sequences` doesn't select file before first click
- Add options for "Show Duplicated" in `Total Sequences`
- Add `exportOverrepresented()`
- **plotSequenceContent**
    - Hide colour information from `hoverinfo`

## Hien

- Genome build information for `gcTheoretical`
- Vignette: 
    - How to create a new object of class `TheoreticalGC` for other organisms
        - Fasta to BSgenome (just point to BSgeneome help)
        - Generate GC content using Mike Love's scripts/package

## Steve

- **Add methods for zero results for all plots, based on plotNContent()**
- Validation Functions for `fastqcData` and `fastqcDataLists` objects
- Build unit tests for classes and methods
- **Test on outlier FastQC reports (1 sequence etc)**
- Finish Vignette
    - Add summary tables
- `Per_base-sequence_quality` colour scheme is incorrect when all > 30
- **plotOverrepresentedSummary**
    - Check plot for no overrepresented sequences: *Temporary version added...*    
    - Include in default report
- Remove RProj files etc
- Add inst/NEWS

# Important But Not Pressing

- Test the output of `runFastQC` now it has been changed to a `FastqcFileList`
- Set character methods for `runFastQC` & redefine this function as a method for class `FastqFileList`

# Future Plans

- Include the capacity for selective overwriting in `runFastQC`
- Calculate Ranks on each module and summarise
- Add a function to merge kMers
- Add Fastq Illumina Filter plot/status
