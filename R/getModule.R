#' @title Retrieve a given module from a Fastqc* Object
#'
#' @description Retrieve a specific module from a Fastqc* object as a
#' data.frame
#'
#' @details
#' This function will return a given module from a Fastqc* object as a
#' data.frame. Note that each module will be it's own unique structure,
#' although all will return a data.frame
#'
#' @param object Can be a \code{FastqcFile}, \code{FastqcFileList},
#' \code{FastqcData}, \code{fastqcDataList}, or simply a \code{character} vector
#' of paths to fastqc files
#' @param module The requested module as contained in a FastQC report. Possible
#' values are \code{Summary}, \code{Basic_Statistics},
#' \code{Per_base_sequence_quality}, \code{Per_tile_sequence_quality},
#' \code{Per_sequence_quality_scores}, \code{Per_base_sequence_content},
#' \code{Per_sequence_GC_content}, \code{Per_base_N_content},
#' \code{Sequence_Length_Distribution}, \code{Sequence_Duplication_Levels},
#' \code{Overrepresented_sequences}, \code{Adapter_Content},
#' \code{Kmer_Content}, \code{Total_Deduplicated_Percentage}.
#' Note that spelling and capitalisation is exactly as contained within a
#' FastQC report, with the exception that spaces have been converted to
#' underscores. Partial matching is implemented for this argument.
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#' @include FastqcFile.R
#' @include FastqcFileList.R
#' @include FastqcDataList.R
#'
#' @return A single \code{tibble} containing module-level information
#' from all FastQC reports contained in the Fastqc* object.
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Extract the Summary module, which corresponds to the PASS/WARN/FAIL flags
#' getModule(fdl, "Summary")
#'
#' # The Basic_Statistics module corresponds to the table at the top of each
#' # FastQC report
#' getModule(fdl, "Basic_Statistics")
#'
#' @docType methods
#'
#' @import tibble
#'
#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastqcData", function(object, module){

    ## This is the list of defined modules in the object specification
    allMods <-  c(
        "Summary", "Basic_Statistics", "Per_base_sequence_quality",
        "Per_tile_sequence_quality", "Per_sequence_quality_scores",
        "Per_base_sequence_content", "Per_sequence_GC_content",
        "Per_base_N_content", "Sequence_Length_Distribution",
        "Sequence_Duplication_Levels", "Overrepresented_sequences",
        "Adapter_Content", "Kmer_Content", "Total_Deduplicated_Percentage"
    )

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)

    ## Get the Filename for the underlying Fastq file
    nm <- fileName(object)

    ## Extract the actual module
    df <- slot(object, module)

    ## The Total_Deduplicated_Percentage needs to be handled separately
    ## as this returns a double which needs to be coerced into a df
    if (module == "Total_Deduplicated_Percentage") {
        if (is.na(df)) {
            df <- tibble()
        }
        else {
            df <- tibble(Total_Deduplicated_Percentage = df)
        }
    }

    ## If we have an empty data.frame, just return that here with a message
    if (nrow(df) == 0) {
        message(module, " missing from ", nm)
        return(df)
    }

    ## Make sure the first column is the name of the underlying fastq file
    cols <- colnames(df)
    df$Filename <- nm

    ## Return the tibble / data.frame
    df[unique(c("Filename", cols))]

})

#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastqcDataList", function(object, module){

    ## This is the list of defined modules on the object specification
    allMods <-  c("Summary", "Basic_Statistics", "Per_base_sequence_quality",
                  "Per_tile_sequence_quality", "Per_sequence_quality_scores",
                  "Per_base_sequence_content", "Per_sequence_GC_content",
                  "Per_base_N_content", "Sequence_Length_Distribution",
                  "Sequence_Duplication_Levels", "Overrepresented_sequences",
                  "Adapter_Content", "Kmer_Content",
                  "Total_Deduplicated_Percentage")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)

    ## Extract module as a list
    allDfs <- lapply(object, getModule, module = module)

    ## Now handle each module as appropriate. All can simply be joined using
    ## dplyr::bind_rows()
    dplyr::bind_rows(allDfs)

})

#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastqcFile", function(object, module){

    ## This is the list of defined modules on the object specification
    allMods <-  c("Summary", "Basic_Statistics", "Per_base_sequence_quality",
                  "Per_tile_sequence_quality", "Per_sequence_quality_scores",
                  "Per_base_sequence_content", "Per_sequence_GC_content",
                  "Per_base_N_content", "Sequence_Length_Distribution",
                  "Sequence_Duplication_Levels", "Overrepresented_sequences",
                  "Adapter_Content", "Kmer_Content",
                  "Total_Deduplicated_Percentage")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)
    ## Obtain a FastqcData object then extract the module
    fqcData <- getFastqcData(object)
    getModule(fqcData, module)

})

#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastqcFileList", function(object, module){

    ## This is the list of defined modules on the object specification
    allMods <-  c("Summary", "Basic_Statistics", "Per_base_sequence_quality",
                  "Per_tile_sequence_quality", "Per_sequence_quality_scores",
                  "Per_base_sequence_content", "Per_sequence_GC_content",
                  "Per_base_N_content", "Sequence_Length_Distribution",
                  "Sequence_Duplication_Levels", "Overrepresented_sequences",
                  "Adapter_Content", "Kmer_Content",
                  "Total_Deduplicated_Percentage")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)
    ## Obtain a FastqcDataList object then extract the module
    fqcData <- getFastqcData(object)
    getModule(fqcData, module)

})

#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "character", function(object, module){

    ## This is the list of defined modules on the object specification
    allMods <-  c("Summary", "Basic_Statistics", "Per_base_sequence_quality",
                  "Per_tile_sequence_quality", "Per_sequence_quality_scores",
                  "Per_base_sequence_content", "Per_sequence_GC_content",
                  "Per_base_N_content", "Sequence_Length_Distribution",
                  "Sequence_Duplication_Levels", "Overrepresented_sequences",
                  "Adapter_Content", "Kmer_Content",
                  "Total_Deduplicated_Percentage")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)
    ## Obtain a FastqcData/List object then extract the module
    fqcData <- getFastqcData(object)
    getModule(fqcData, module)

})
