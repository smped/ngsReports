#' @title Get the Adapter Content information
#'
#' @description Retrieve the Adapter Content module from one or more FastQC
#' reports
#'
#' @details It is currently assumed that all FastQC reports contain information
#' about the same adapters. If this was varied during the preparation of the
#' reports, this function will fail.
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList},
#' \code{FastqcData}, \code{fastqcDataList}, or simply a \code{character}
#' vector of paths to fastqc files.
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#' @include FastqcFile.R
#' @include FastqcFileList.R
#' @include FastqcDataList.R
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
#' # Return a tibble/tibble with the raw information
#' Adapter_Content(fdl)
#'
#' @import tibble
#'
#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod(
    f = "Adapter_Content",
    signature = "FastqcDataList",
    definition = function(object){
        df <- lapply(object@.Data, Adapter_Content)
        nulls <- vapply(df,
                        function(x){
                            length(x) == 0
                        }, logical(1))
        if (sum(nulls) > 0) message(
            sprintf("The Adapter_Content module was missing from:\n%s",
                    paste(path(object)[nulls], sep = "\n"))
        )
        dplyr::bind_rows(df)
    })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod(
    f = "Adapter_Content",
    signature = "FastqcData",
    definition = function(object){
        df <- object@Adapter_Content
        if (length(df)) { # Check there is data in the module
            ## Add a Filename column if there is any data
            df$Filename <- fileName(object)
            dplyr::select(df, "Filename", tidyselect::everything())
        }
        else {# Otherwise return the blank data.frame
            df
        }
    })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod(
    f = "Adapter_Content",
    signature = "FastqcFile",
    definition = function(object){
        object <- getFastqcData(object)
        Adapter_Content(object)
    })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod(
    f = "Adapter_Content",
    signature = "FastqcFileList",
    definition = function(object){
        object <- getFastqcData(object)
        Adapter_Content(object)
    })

#' @export
#' @rdname Adapter_Content
#' @aliases Adapter_Content
setMethod(
    f = "Adapter_Content",
    signature = "character",
    definition = function(object){
        object <- getFastqcData(object)
        Adapter_Content(object)
    })
