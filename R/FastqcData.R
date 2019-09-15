#' @title The FastqcData Object Class
#'
#' @description The FastqcData Object Class
#'
#' @details This object class is the main object required for generating plots
#' and tables. Instantiation will first test for a compressed file (or
#' extracted directory) with the correct data structure, and will then parse
#' all the data into R as a \code{FastqcData} object. FastQC modules are
#' contained as individual slots, which can be viewed using \code{slotNames}.
#'
#' Individual modules can be returned using the function \code{getModule()}
#' and specifying which module is required. See \code{\link{getModule}} for
#' more details.
#'
#' @slot Summary Summary of PASS/WARN/FAIL status for each module
#' @slot Basic_Statistics The Basic_Statistics table from the top of a FastQC
#' html report
#' @slot Per_base_sequence_quality The underlying data from the
#' Per_base_sequence_quality module
#' @slot Per_sequence_quality_scores The underlying data from the
#' Per_sequence_quality_scores module
#' @slot Per_base_sequence_content The underlying data from the
#' Per_base_sequence_content module
#' @slot Per_sequence_GC_content The underlying data from the
#' Per_sequence_GC_content module
#' @slot Per_base_N_content The underlying data from the
#' Per_base_N_content module
#' @slot Sequence_Length_Distribution The underlying data from the
#' Sequence_Length_Distribution module
#' @slot Sequence_Duplication_Levels The underlying data from the
#' Sequence_Duplication_Levels module
#' @slot Overrepresented_sequences The underlying data from the
#' Overrepresented_sequences module
#' @slot Adapter_Content The underlying data from the Adapter_Content module
#' @slot Kmer_Content The underlying data from the Kmer_Content module
#' @slot Total_Deduplicated_Percentage Estimate taken from the plot data for
#' Sequence_Duplication_Levels. Only included in later versions of FastQC
#' @slot version The version of FastQC used for generation of the report (if
#' available)
#' @slot path Path to the FastQC report#'
#'
#' @return An object of class FastqcData
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)[1]
#'
#' # Load the FASTQC data as a FastqcData object
#' fd <- FastqcData(fl)
#' fd
#'
#' @include validationFunctions.R
#'
#' @rdname FastqcData
#' @aliases FastqcData-class
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
        version = "character",
        path = "character"
    )
)
setValidity("FastqcData", .isValidFastqcData)

#' @param x Path to a single zip archive or extracted folder for a individual
#' FastQC report.
#' @rdname FastqcData
#' @export
FastqcData <- function(x){
    stopifnot(!is.null(x))
    fl <- .FastqcFile(x)
    as(fl, "FastqcData")
}

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
