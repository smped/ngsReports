#' @title Get the summary information from Fastqc Files
#'
#' @description Read the information from the \code{summary.txt} files in each
#' .FastqcFile
#'
#' @details
#' This simply extracts the summary of PASS/WARN/FAIL status for every module
#' as defined by the tool FastQC for each supplied file.
#'
#' @param object Can be a \code{FastqcData}, \code{FastqcDataList} object or
#' a vector of paths to unparsed FastQC reports.
#'
#' @return A \code{tibble} containing the PASS/WARN/FAIL status for each
#' module, as defined in a FastQC report.
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Return a tibble/tibble with the raw information
#' getSummary(fdl)
#'
#' @importFrom utils unzip
#' @import tibble
#'
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", ".FastqcFile", function(object){
    modules <- c(
        "Basic Statistics",
        "Per base sequence quality",
        "Per tile sequence quality",
        "Per sequence quality scores",
        "Per base sequence content",
        "Per sequence GC content",
        "Per base N content",
        "Sequence Length Distribution",
        "Sequence Duplication Levels",
        "Overrepresented sequences",
        "Adapter Content",
        "Kmer Content"
    )
    path <- path(object)
    if (isCompressed(path, type = "zip")) {
        ##Get the internal path within the zip archive
        if (!file.exists(path)) stop("The zip archive can not be found.")
        fl <- file.path(gsub(".zip$", "", basename(path)), "summary.txt")
        ## Check the required file exists
        allFiles <- unzip(path, list = TRUE)$Name
        if (!fl %in% allFiles)
            stop("summary.txt is missing from the zip archive")
        ## Open the connection & read all lines
        uz <- unz(path,fl)
        summaryLines <- readLines(uz)
        close(uz)

        ## Form the output
        summaryData <- .splitByTab(summaryLines, firstRowToNames = FALSE)
    }
    else{
        ## The existence of this file will have been checked at object
        ## instantion
        ## Check in case it has been deleted post-instantiation though
        fl <- file.path(path, "summary.txt")
        if (!file.exists(fl)) stop("'summary.txt' could not be found.")
        summaryData <- suppressMessages(
            readr::read_tsv(fl, col_names = FALSE)
        )
    }
    colnames(summaryData) <- c("Status", "Category", "Filename")
    if (!any(modules %in% summaryData$Category))
        stop("summary.txt contained none of the expected modules.")
    ## Return output as a tibble
    as_tibble(summaryData)
})

#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "ANY", function(object){
    ## This will error if the file doesn't exist, or if it's not in the
    ## correct format (i.e a FastQC output)
    object <- FastqcDataList(object)
    if (length(object) == 1) object <- object[[1]]
    getSummary(object)
})

#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "FastqcData", function(object){object@Summary})

#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "FastqcDataList", function(object){
    df <- lapply(object@.Data, getSummary)
    dplyr::bind_rows(df)
})
