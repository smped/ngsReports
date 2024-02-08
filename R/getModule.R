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
#' @param object Can be a `FastqcData`, `fastqcDataList`, or simply
#' a `character` vector of paths
#' @param module The requested module as contained in a FastQC report. Possible
#' values are `Summary`, `Basic_Statistics`,
#' `Per_base_sequence_quality`, `Per_tile_sequence_quality`,
#' `Per_sequence_quality_scores`, `Per_base_sequence_content`,
#' `Per_sequence_GC_content`, `Per_base_N_content`,
#' `Sequence_Length_Distribution`, `Sequence_Duplication_Levels`,
#' `Overrepresented_sequences`, `Adapter_Content`,
#' `Kmer_Content`, `Total_Deduplicated_Percentage`.
#' Note that spelling and capitalisation is exactly as contained within a
#' FastQC report, with the exception that spaces have been converted to
#' underscores. Partial matching is implemented for this argument.
#'
#' @include FastqcData.R
#' @include AllGenerics.R
#' @include FastqcFile.R
#' @include FastqcDataList.R
#'
#' @return A single `tibble` containing module-level information
#' from all FastQC reports contained in the Fastqc* object.
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
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
#' @importFrom methods slotNames slot
#'
#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastqcData", function(object, module){

    ## This is the list of defined modules from the object specification
    allMods <-  slotNames("FastqcData")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)

    ## Get the Filename for the underlying Fastq file
    nm <- fqName(object)

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

    ## This is the list of defined modules from the object specification
    allMods <-  slotNames("FastqcData")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)

    ## Extract module as a list
    allDfs <- lapply(object, getModule, module = module)

    ## Now handle each module as appropriate. All can simply be joined using
    ## dplyr::bind_rows()
    allDfs <- dplyr::bind_rows(allDfs)
    as_tibble(allDfs)

})

#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "ANY", function(object, module){

    ## Obtain a FastqcData/List object then extract the module
    ## This will fail if an invalid file path is specified
    fqcData <- FastqcDataList(object)
    if (length(fqcData) == 1) fqcData <- fqcData[[1]]
    getModule(fqcData, module)

})

#' @importFrom tibble as_tibble_row
#' @importFrom dplyr bind_cols
#' @importFrom methods slotNames slot
#'
#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastpData", function(object, module){

    ## This is the list of defined modules from the object specification
    allMods <-  slotNames("FastpData")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)

    ## Get the Filename for the underlying Fastq file
    isPaired <- object@paired
    nm <- c(Filename = basename(path(object)), fqName(object))
    nm_tbl <- as_tibble_row(nm)

    ## Extract the actual module
    mod <- slot(object, module)

    if (module == "Summary") {
        mod <- lapply(mod, function(x) bind_cols(nm_tbl, x))
    }
    if (module %in% c("Before_filtering", "After_filtering")) {
        mod$read1 <- bind_cols(nm_tbl[c("Filename", "read1")], mod$read1)
        names(mod$read1)[names(mod$read1) == "read1"] <- "fqName"
        if (isPaired) {
            mod$read2 <- bind_cols(nm_tbl[c("Filename", "read2")], mod$read2)
            names(mod$read2)[names(mod$read2) == "read2"] <- "fqName"
        }
    }
    if (module %in% c("Adapters", "Duplication", "Insert_size")) {
        mod <- bind_cols(nm_tbl, mod)
    }

    mod

})

#' @export
#' @rdname getModule
#' @aliases getModule
setMethod("getModule", "FastpDataList", function(object, module){

    ## This is the list of defined modules from the object specification
    allMods <-  slotNames("FastpData")

    ## Make sure we have asked for a valid module
    module <- match.arg(module, allMods)

    ## Extract module as a list
    anyPaired <- any(vapply(object, function(x) x@paired, logical(1)))
    allMods <- lapply(object, getModule, module = module)
    out <- list()

    if (module == "Summary") {
        modNames <- c("Before_filtering", "After_filtering", "Filtering_result")
        out <- lapply(
            modNames,
            function(mod) bind_rows(lapply(allMods, function(x) x[[mod]]))
        )
        names(out) <- modNames
    }
    if (module %in% c("Before_filtering", "After_filtering")) {
        r1 <- bind_rows(lapply(allMods, function(x) x[["read1"]]))
        out[["read1"]] <- r1
        if (anyPaired) {
            r2 <- bind_rows(lapply(allMods, function(x) x$read2))
            out[["read2"]] <- r2
        }
    }
    if (module %in% c("Adapters", "Duplication", "Insert_size")) {
        out <- bind_rows(allMods)
    }
    if (module %in% c("paired", "command", "version", "path")) {
        vals <- unlist(allMods)
        names(vals) <- c()
        out <- tibble(Filename = basename(path(object)))
        out[[module]] <- vals
    }
    out

})
