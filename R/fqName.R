#' @title Return the Underlying Fastq File Names from Fastqc/Fastp Objects
#'
#' @description Return the Underlying Fastq File Names from Fastqc/Fastp Objects
#'
#' @param object An object able to extract an Fastq name from
#' @param value Replacement value for fqName
#'
#' @return Returns the names of the Fastq files the FastQC report was
#' generated from, without any preceding directories.
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#' fqName(fdl)
#'
#' nm <- paste0(letters[seq_along(fdl)], ".fq")
#' fqName(fdl) <- nm
#' fqName(fdl)
#'
#'
#' @export
#' @rdname fqName-methods
setGeneric("fqName", function(object){standardGeneric("fqName")})
#' @export
#' @rdname fqName-methods
setMethod("fqName", "ANY", function(object){
    cl <- class(object)
    message("Method 'fqName' not implemented for objects of class ", cl)
})
#' @export
#' @name fqName
#' @aliases fqName,FastqcData-method
#' @rdname fqName-methods
setMethod("fqName", "FastqcData", function(object){
    object@Summary$Filename[1]
})
#' @export
#' @name fqName
#' @aliases fqName,FastqcDataList-method
#' @rdname fqName-methods
setMethod("fqName", "FastqcDataList", function(object){
    vapply(object@.Data, fqName, character(1))
})
#' @export
#' @rdname fqName-methods
setGeneric(
    "fqName<-", signature = "object",
    function(object, value) standardGeneric("fqName<-")
)
#' @export
#' @name fqName<-
#' @aliases fqName<-,FastqcData-method
#' @rdname fqName-methods
setReplaceMethod("fqName", signature = "FastqcData", function(object, value){
    stopifnot(length(value) == 1 & is.character(value))
    n <- nrow(object@Summary)
    object@Summary$Filename <- rep(value, n)
    object@Basic_Statistics$Filename <- value
    object
})
#' @export
#' @name fqName<-
#' @aliases fqName<-,FastqcDataList-method
#' @rdname fqName-methods
setReplaceMethod("fqName", signature = "FastqcDataList", function(object, value){
    stopifnot(length(value) == length(object))
    stopifnot(is.character(value))
    out <- lapply(
        seq_along(value),
        function(i, x = object){
            fqName(x[[i]]) <- value[i]
            x[[i]]
        }
    )
    as(out, "FastqcDataList")
})
#' @export
#' @name fqName
#' @aliases fqName,FastpData-method
#' @rdname fqName-methods
setMethod("fqName", "FastpData", function(object){
    isPaired <- object@paired
    cmd <- object@command
    r1 <- basename(gsub(".+(-i|--in1) *([^ ]+) .+", "\\2", cmd))
    out <- c(read1 = r1)
    if (isPaired) {
        r2 <- basename(gsub(".+(-I|--in2) *([^ ]+) .+", "\\2", cmd))
        out <- c(out, read2 = r2)
    }
    ## If names couldn't be extracted from the command, set as NA
    any_na <- out == cmd
    out[any_na] <- NA_character_
    out
})
#' @export
#' @name fqName
#' @aliases fqName,FastpDataList-method
#' @rdname fqName-methods
setMethod("fqName", "FastpDataList", function(object){
    allPaired <- all(vapply(object, function(x) x@paired, logical(1)))
    if (allPaired) {
        fqName <- t(vapply(object@.Data, fqName, character(2)))
    } else {
        fqName <- vapply(object@.Data, function(x) fqName(x)[1], character(1))
    }
    fqName
})
