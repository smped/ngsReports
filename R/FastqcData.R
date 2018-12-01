#' @title The FastqcData Object Class
#'
#' @description The FastqcData Object Class
#'
#' @return An object of class FastqcData
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcData object
#' # As this is the underlying structure for a FastqcDataList, an object of
#' # FastqcData will only be returned from an individual file.
#' fd <- getFastqcData(fileList[1])
#' fd
#'
#' @include validationFunctions.R
#'
#' @slot ... this can either be a single character vector of paths to FASTQC
#' files, or several instances of FastqcFile objects
setClass(
    "FastqcData",
    slots = c(
        Summary = "data.frame",
        Basic_Statistics = "data.frame",
        Per_base_sequence_quality = "data.frame",
        Per_tile_sequence_quality = "data.frame",
        Per_sequence_quality_scores = "data.frame",
        Per_base_sequence_content = "data.frame",
        Per_sequence_GC_content = "data.frame",
        Per_base_N_content = "data.frame",
        Sequence_Length_Distribution = "data.frame",
        Sequence_Duplication_Levels = "data.frame",
        Overrepresented_sequences = "data.frame",
        Adapter_Content = "data.frame",
        Kmer_Content = "data.frame",
        Total_Deduplicated_Percentage = "numeric",
        Version = "character",
        path = "character"
    )
)
setValidity("FastqcData", .isValidFastqcData)

## The show method doesn't need exporting
setMethod(
    "show",
    "FastqcData",
    function(object){
        cat("FastqcData for", object@Basic_Statistics$Filename, "\n")
        cat("Source Fastq file contains",
            scales::comma(object@Basic_Statistics$Total_Sequences),
            "reads.\n")
        cat("Source FastQC file is located in", object@path)
    }
)
