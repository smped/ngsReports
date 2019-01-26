#' @title Get the Sequence Duplication Levels information
#'
#' @description Retrieve the Sequence Duplication Levels module from one or
#' more FastQC reports
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList},
#' \code{FastqcData}, \code{fastqcDataList}, or simply a \code{character}
#' vector of paths to fastqc files
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#'
#' @return A single \code{tibble} containing all information combined from all
#' supplied FastQC reports
#'
#'@examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Print the Sequence Duplication Levels
#' Sequence_Duplication_Levels(fdl)
#'
#' @docType methods
#'
#'
#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcData", function(object){
    df <- object@Sequence_Duplication_Levels
    if (length(df)) {# Check there is data in the module
        ## Add a Filename column if there is any data
        df$Filename <- fileName(object)
        dplyr::select(df, "Filename", tidyselect::everything())
    }
    else {# Otherwise return the blank data.frame
        df
    }
}
)
#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcDataList", function(object){
    df <- lapply(object@.Data, Sequence_Duplication_Levels)
    # Find any for which no data was returned
    nulls <- vapply(df, function(x){length(x) == 0}, logical(1))
    if (sum(nulls) > 0)
        message(
            sprintf(
                "Sequence_Duplication_Levels module missing from %s\n",
                paste(path(object)[nulls], sep = "\n"))
        )
    dplyr::bind_rows(df)
}
)
#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcFile", function(object){
    object <- getFastqcData(object)
    Sequence_Duplication_Levels(object)
}
)
#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "FastqcFileList", function(object){
    object <- getFastqcData(object)
    Sequence_Duplication_Levels(object)
}
)
#' @export
#' @rdname Sequence_Duplication_Levels
#' @aliases Sequence_Duplication_Levels
setMethod("Sequence_Duplication_Levels", "character", function(object){
    object <- getFastqcData(object)
    Sequence_Duplication_Levels(object)
}
)
