#' @importFrom utils unzip
#'
.isValidFastqcFile <- function(object){

    path <- path(object)
    if (length(path) != 1) {
        warning("Only a single path can be specified")
        return(FALSE)
    }

    comp <- isCompressed(path, type = "zip") # This includes a file.exists step
    if (comp) {
        ## List the files
        subFiles <- basename(unzip(path, list = TRUE)$Name)
    }
    else{
        ## Check the file is a directory
        chk <- checkmate::testDirectoryExists(path)
        if (!chk) {
            warning("The supplied file is not a directory")
            return(FALSE)
        }
        subFiles <- list.files(path)
    }
    if (any(!c("fastqc_data.txt", "summary.txt") %in% subFiles)) {
        warning("The required files are missing from the supplied path/file")
        return(FALSE)
    }
    ## If it checks out up to here, we're good to go
    TRUE
}

#' @importFrom methods slotNames
.isValidFastqcData <- function(object){
    ## At minimum, this should contain slots for Summary, Basic_Statistics,
    ## version & path. The remainder of the modules may be missing, although
    ## they will stil be present, just as empty data.frame objects.
    reqSlots <- c("Summary", "Basic_Statistics", "version", "path")
    all(reqSlots %in% slotNames(object))
}

#'
.isValidFastpFile <- function(object){

    path <- path(object)
    if (length(path) != 1) {
        warning("Only a single path can be specified")
        return(FALSE)
    }
    expected <- c("{", "\t\"summary\": {")
    head <- readLines(path, 2)
    all(expected == head)
}

#' @importFrom methods slotNames
.isValidFastpData <- function(object){
    reqSlots <- c(
        "Summary", "Adapters", "Duplication", "Insert_size",
        "Before_filtering", "After_filtering"
    )
    all(reqSlots %in% slotNames(object))
}

#' @importFrom methods is
.isValidFastqcDataList <- function(object){
    ## This is very rudimentary & may need more thought
    cls <- vapply(object, is, logical(1), "FastqcData")
    all(cls)
}

#' @importFrom methods is
.isValidFastpDataList <- function(object){
    ## This is very rudimentary & may need more thought
    cls <- vapply(object, is, logical(1), "FastpData")
    all(cls)
}

.isValidPwf <- function(object){

    vals <- getColours(object)
    if (length(vals) != 4) return(FALSE) # Ensure all have been set
    if (any(substr(vals, 1, 1) != "#")) return(FALSE) # Start with #
    if (any(!nchar(vals) %in% c(7, 9))) return(FALSE) # RGB length 7 or 9
    if (any(grepl("[G-Z]", vals))) return(FALSE) # Invald Hex values

    TRUE

}

#' @importFrom methods slot
.isValidTheoreticalGC <- function(object){

    ## Check the colnames of the metaData
    reqCols <- c("Name", "Genome", "Transcriptome")
    stopifnot(all(reqCols %in% colnames(object@mData)))
    gn <- object@mData$Genome
    stopifnot(is.logical(gn))
    tr <- object@mData$Transcriptome
    stopifnot(is.logical(tr))

    ## Check Genome & Transcriptomes called as TRUE match the metadata exactly
    if (any(gn)) {
        stopifnot(all(object@mData$Name[gn] %in% colnames(object@Genome)))
        stopifnot(all(colnames(object@Genome)[-1] %in% object@mData$Name[gn]))
        ## Check all content is given as a frequency
        gFreqs <- vapply(
            object@Genome[-1], function(x){sum(x > 1 | x < 0)}, integer(1)
        )
        stopifnot(all(gFreqs == 0))
    }
    if (any(tr)) {
        stopifnot(
            all(object@mData$Name[tr] %in% colnames(object@Transcriptome))
        )
        stopifnot(
            all(colnames(object@Transcriptome)[-1] %in% object@mData$Name[tr])
        )
        ## Check the frequencies again
        tFreqs <- vapply(
            object@Transcriptome[-1],
            function(x){sum(x > 1 | x < 0)},
            integer(1)
        )
        stopifnot(all(tFreqs == 0))
    }

    TRUE

}
