#' @title Get the Kmer Content information
#'
#' @description Retrieve the Kmer Content module from one or more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList},
#' \code{FastqcData}, \code{fastqcDataList}, or a \code{character} vector
#' of file paths
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{tibble} containing all information combined from all
#' supplied FastQC reports
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Print the Kmer_Content
#' Kmer_Content(fdl)
#'
#' @docType methods
#'
#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcData", function(object){
    df <- dplyr::mutate(object@Kmer_Content)
    if (length(df)) {# Check there is data in the module
        ## Add a Filename column if there is any data
        df$Filename <- fileName(object)
        dplyr::select(df, "Filename", tidyselect::everything())
    }
    else {# Otherwise return the blank data.frame
        df
    }
})

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcDataList", function(object){
    df <- lapply(object@.Data, Kmer_Content)
    nulls <- vapply(df,
                    function(x){
                        length(x) == 0
                    }, logical(1))
    if (sum(nulls) > 0) message(
        sprintf(
            "The Kmer_Content module was missing from %s\n",
            paste(path(object)[nulls], sep = "\n")
        )
    )
    dplyr::bind_rows(df)
})

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcFile", function(object){
    object <- getFastqcData(object)
    Kmer_Content(object)
})

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "FastqcFileList", function(object){
    object <- getFastqcData(object)
    Kmer_Content(object)
})

#' @export
#' @rdname Kmer_Content
#' @aliases Kmer_Content
setMethod("Kmer_Content", "character", function(object){
    object <- getFastqcData(object)
    Kmer_Content(object)
})
