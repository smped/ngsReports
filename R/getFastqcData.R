#' Get all data from FastQC files
#'
#' @description Read the information from the \code{fastqc_data.txt} files in
#' each FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList, or paths to files
#'
#' @return An object of \code{FastqcData} or a \code{FastqcDataList}
#'
#' @include AllGenerics.R
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' @export
#' @rdname getFastqcData-methods
setGeneric("getFastqcData", function(object){standardGeneric("getFastqcData")})

#' @importFrom utils unzip
#' @name getFastqcData
#' @aliases getFastqcData,FastqcFile-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "FastqcFile", function(object){

    ## Setup the path to the file & use tryCatch to ensure it exists
    ## when testing to see if it is compressed
    path <- path(object)
    comp <- tryCatch(isCompressed(path, type = "zip"))

    if (comp) {

        ## Get the internal path within the zip archive
        fl <- file.path(gsub(".zip$", "", basename(path)), "fastqc_data.txt")

        ## Check the required file exists within the file
        allFiles <- unzip(path, list = TRUE)$Name
        stopifnot(fl %in% allFiles)

        ## # Open the connection & read the 12 lines
        uz <- unz(path, fl)
        fqcLines <- readLines(uz)
        close(uz)

    }
    else{

        ## The existence of this file will have been checked at
        ## instantiation of the FastqcFile. Check in case it has been
        ## deleted post-instantiation
        fl <- file.path(path, "fastqc_data.txt")
        if (!file.exists(fl)) stop("'fastqc_data.txt' could not be found.")
        fqcLines <- readLines(fl)

    }

    ## Modified from Repitools::readFastQC
    ## Remove any '#' symbols
    fqcLines <- gsub("#", "", fqcLines)
    ## Remove any lines which specify '>>END_MODULE'
    fqcLines <- fqcLines[!grepl(">>END_MODULE", fqcLines)]

    ## The FastQC version NUMBER
    vers <- NA_character_ # Default to missing
    if (grepl("[Ff][Aa][Ss][Tt][Qq][Cc]\\t", fqcLines[1])) {
        vers <- gsub("[Ff][Aa][Ss][Tt][Qq][Cc]\\t(.+)", "\\1", fqcLines[1])
        ## Remove the version line for easier module identification
        fqcLines <- fqcLines[-1]
    }

    ## Setup the module names
    modules <- grep("^>>", fqcLines, value = TRUE)
    ## Extract the text after the '>>' & before the tab
    modules <- gsub("^>>(.+)\\t.+", "\\1", modules)
    modules <- gsub(" ", "_", modules) # Add underscores

    ## Define the standard modules
    reqModules <- c(
        "Basic_Statistics",
        "Per_base_sequence_quality",
        "Per_tile_sequence_quality",
        "Per_sequence_quality_scores",
        "Per_base_sequence_content",
        "Per_sequence_GC_content",
        "Per_base_N_content",
        "Sequence_Length_Distribution",
        "Sequence_Duplication_Levels",
        "Overrepresented_sequences",
        "Adapter_Content",
        "Kmer_Content"
    )
    ## Check that at least one of the standard modules is present
    if (!any(modules %in% reqModules))
        stop("None of the standard modules were found in the data.")

    ## Split the data based on the '>>' pattern, which indicates the
    ## beginning of a new module
    fqcLines <- split(fqcLines, cumsum(grepl("^>>", fqcLines)))

    ## Assign the module names
    names(fqcLines) <- vapply(fqcLines, function(x){
        ## Get the first line before the tab separator
        nm <- gsub(">>(.+)\\t.+", "\\1", x[[1]])
        ## Replace the space with underscore
        gsub(" ", "_", nm)
    }, character(1))
    ## Remove the module name (i.e. the 1st value) from each module's data
    fqcLines <- lapply(fqcLines, function(x){x[-1]})

    ## Define the output to have the same structure as fastqcData,
    ## with an additional slot to account for the Sequence Duplication
    ## Levels output. Initialise an empty list based on the standard
    ## modules after checking for the Total_Deduplicated Percentage.
    allModules <- c(reqModules, modules, "Total_Deduplicated_Percentage")
    allModules <- unique(allModules)
    out <- vector("list", length(allModules))
    names(out) <- allModules

    ## Get the Modules
    out[["Basic_Statistics"]] <- .getBasicStatistics(fqcLines)
    out[["Per_base_sequence_quality"]] <- .getPerBaseSeqQuals(fqcLines)
    out[["Per_tile_sequence_quality"]] <- .getPerTileSeqQuals(fqcLines)
    out[["Per_sequence_quality_scores"]] <- .getPerSeqQualScores(fqcLines)
    out[["Per_base_sequence_content"]] <- .getPerBaseSeqCont(fqcLines)
    out[["Per_sequence_GC_content"]] <- .getPerSeqGcCont(fqcLines)
    out[["Per_base_N_content"]] <-  .getPerBaseNCont(fqcLines)
    out[["Sequence_Length_Distribution"]] <- .getSeqLengthDist(fqcLines)
    out[["Overrepresented_sequences"]] <- .getOverrepSeq(fqcLines)
    out[["Adapter_Content"]] <- .getAdapterCont(fqcLines)
    out[["Kmer_Content"]] <- .getKmerCont(fqcLines)

    ## Get the Sequence Duplication Levels
    Sequence_Duplication_Levels <- .getSeqDuplicationLevels(fqcLines)
    dupMods <- names(Sequence_Duplication_Levels)
    out[dupMods] <- Sequence_Duplication_Levels[dupMods]

    ##Get the summary
    Summary <- getSummary(object)

    args <- c(
        list(
            Class = "FastqcData",
            path = path(object),
            Version = vers,
            Summary = Summary),
        out
    )

    do.call("new", args)

})

#' @name getFastqcData
#' @aliases getFastqcData,NULL-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "NULL", function(object){
    if (is.null(object)) stop("No files have been provided.")
}
)
#' @name getFastqcData
#' @aliases getFastqcData,character-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "character", function(object){
    x <- FastqcFileList(object)
    n <- length(x)
    if (n == 1) x[[1]]
    getFastqcData(x)
}
)
#' @name getFastqcData
#' @aliases getFastqcData,FastqcFileList-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "FastqcFileList", function(object){
    fqc <- lapply(object@.Data, getFastqcData)
    new("FastqcDataList", fqc)
}
)

