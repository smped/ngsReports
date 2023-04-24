#' @title The FastpData Object Class
#'
#' @description The FastpData Object Class
#' `r lifecycle::badge("experimental")`
#'
#' @details This object class is the main object required for generating plots
#' and tables. Instantiation will first check for a .json file with the correct
#' data structure, and will then parse all the data into R as a `FastpData`
#' object. Fastp modules are contained as individual slots, which can be viewed
#' using `slotNames`. Sub-modules are also contained within many larger modules
#' with modules being based on the sections within a fastp html report
#'
#' Individual modules can be returned using the function `getModule()`
#' and specifying which module/sub-module is required. See [getModule()] for
#' more details.
#'
#' @slot Summary Contains three submodules 1) Before_filtering, 2)
#' After_filtering and 3) Filtering_result. All values presented in the initial
#' table for individual fastp reports are contained in other sections of the
#' report
#' @slot Adapters Contains a tibble with all data from this module
#' @slot Duplication Contains a tibble with all duplication results
#' @slot Insert_size Contains a tibble with all insert size estimates
#' @slot Before_filtering,After_filtering The modules can be selected for either
#' Read1 or Read2
#' @slot paired logical(1) indicating whether the file is from paired-end
#' sequencing
#' @slot command character(1) with the executed command
#' @slot version character(1) with the fastp version being used
#' @slot path Path to the Fastp report
#'
#' @return An object of class FastpData
#'
#' @include validationFunctions.R
#'
#' @rdname FastpData
#' @aliases FastpData-class
setClass(
  "FastpData",
  slots = c(
    Summary = "list",
    Adapters = "data.frame",
    Duplication = "data.frame",
    Insert_size = "data.frame",
    Before_filtering = "list",
    After_filtering = "list",
    paired = "logical",
    command = "character",
    version = "character",
    path = "character"
  )
)
setValidity("FastpData", .isValidFastpData)

#' @param x Path to a single zip archive or extracted folder for a individual
#' fastp report.
#' @rdname FastpData
#' @export
FastpData <- function(x){
  stopifnot(!is.null(x))
  fl <- .FastpFile(x)
  as(fl, "FastpData")
}

## The show method doesn't need exporting
setMethod(
  "show",
  "FastpData",
  function(object){
    isPaired <- object@paired
    reads <- object@Summary$Before_filtering$total_reads
    cat("Source Fastp file is located in", object@path)
    cat(
      "\nSource Fastq file(s) contains",
      ifelse(
        isPaired,
        paste(c(scales::comma(reads/2, 1), "paired-end"), collapse = " "),
        paste(c(scales::comma(reads, 1), "single-end"), collapse = " ")
      ),
      "reads.\n"
    )
  }
)

