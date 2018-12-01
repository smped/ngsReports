#' @importFrom utils unzip
#'
.isValidFastqcFile <- function(object){

    path <- path(object)
    if (length(path) != 1) {
        warning("Only a single path can be specified")
        return(FALSE)
    }

    comp <- isCompressed(path, type = "zip") # This includes a file.exists step
    if (comp){
        ## List the files
        subFiles <- basename(unzip(path, list = TRUE)$Name)
    }
    else{
        ## Check the file is a directory
        chk <- checkmate::testDirectoryExists(path)
        if (!chk){
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

.isValidFastqcFileList <- function(object){
    cls <- vapply(object, class, character(1))
    all(cls == "FastqcFile")
}

.isValidFastqcData <- function(object){
    ## At minimum, this should contain slots for Summary, Basic_Statistics,
    ## Version & path. The remainder of the modules may be missing
    reqSlots <- c("Summary", "Basic_Statistics", "Version", "path")
    all(reqSlots %in% slotNames(object))
}

.isValidFastqcDataList <- function(object){
    ## This is very rudimentary & may need more thought
    cls <- vapply(object, class, character(1))
    all(cls == "FastqcData")
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
    if (!all(reqCols %in% colnames(object@mData))) return(FALSE)
    gn <- object@mData$Genome
    if (!is.logical(gn)) return(FALSE)
    tr <- object@mData$Transcriptome
    if (!is.logical(tr)) return(FALSE)

    ## Check Genome & Transcriptomes called as TRUE match the metadata exactly
    if (!all(object@mData$Name[gn] %in% colnames(object@Genome)))
        return(FALSE)
    if (!all(colnames(object@Genome)[-1] %in% object@mData$Name[gn]))
        return(FALSE)
    if (!all(object@mData$Name[tr] %in% colnames(object@Transcriptome)))
        return(FALSE)
    if (!all(colnames(object@Transcriptome)[-1] %in% object@mData$Name[tr]))
        return(FALSE)

    ## Check all content is given as a frequency
    gFreqs <-
        vapply(object@Genome[-1], function(x){sum(x > 1 | x < 0)}, integer(1))
    if (any(gFreqs != 0)) return(FALSE)
    tFreqs <- vapply(
        object@Transcriptome[-1], function(x){sum(x > 1 | x < 0)}, integer(1)
    )
    if (any(tFreqs != 0)) return(FALSE)

    TRUE

}
