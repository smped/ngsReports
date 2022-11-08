#' @title Write fasta of Over-Represented sequences.
#'
#' @description Output overrepresented sequences to disk in fasta format.
#'
#' @details Fasta will contain `Filename`, `Possible Source`,
#' `Percent of total reads`
#'
#' @param x Can be a `FastqcData` or `FastqcDataList`
#' @param path Path to export the fasta file to. Reverts to a default in the
#' working directory if not supplied
#' @param n The number of sequences to output
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param noAdapters logical. Remove any sequences identified as possible
#' adapters or primers by FastQC
#' @param ... Used to pass any alternative patterns to remove from the end of
#' filenames
#'
#' @return Exports to a fasta file, and returns the fasta information invisibly
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
#' # Export the top10 Overrepresented Sequences as a single fasta file
#' faOut <- file.path(tempdir(), "top10.fa")
#' overRep2Fasta(fdl, path = faOut)
#'
#' @docType methods
#'
#' @name overRep2Fasta
#' @rdname overRep2Fasta-methods
#' @export
setGeneric("overRep2Fasta", function(
    x, path, n = 10, labels, noAdapters = TRUE, ...){
    standardGeneric("overRep2Fasta")
})
#' @rdname overRep2Fasta-methods
#' @export
setMethod("overRep2Fasta", signature = "ANY", function(
    x, ...){
    .errNotImp(x)
}
)
#' @rdname overRep2Fasta-methods
#' @export
setMethod("overRep2Fasta", signature = "FastqcData", function(
    x, path, n = 10, labels, noAdapters = TRUE, ...){

    ## Sort out the labels
    df <- getModule(x, "Overrepresented_sequences")
    labels <- .makeLabels(x, labels, ...)
    labels <- labels[names(labels) %in% df$Filename]

    if (missing(path)) path <-paste(
        fqName(x), "top", n, "overrepresented.fa", sep = "_"
    )

    ## Remove any putative adapter or primer sequences
    if (noAdapters) df <- df[!grepl("(Primer|Adapter)", df$Possible_Source),]
    if (nrow(df) == 0) stop("No overrepresented sequences found.")

    ## Now just extract the top n & export
    df <- dplyr::top_n(df, n, Percentage)
    n <- nrow(df) # In case less than the requested n were present
    hdr <- paste0(
        "> ", labels,
        " Sequence", seq_len(n),
        " Count:", df$Count,
        " Length:", nchar(df$Sequence)
    )
    fasta <- paste(hdr, df$Sequence, sep = "\n")
    writeLines(fasta, path)
    invisible(fasta)
}
)
#' @rdname overRep2Fasta-methods
#' @export
setMethod("overRep2Fasta", signature = "FastqcDataList", function(
    x, path, n = 10, labels, noAdapters = TRUE, ...){

    df <- getModule(x, "Overrepresented_sequences")

    if (missing(path)) path <- paste(
        "all_files", "top", n, "overrepresented.fa", sep = "_"
    )

    ## Remove any putative adapter or primer sequences
    if (noAdapters) df <- df[!grepl("(Primer|Adapter)", df$Possible_Source),]
    if (nrow(df) == 0) stop("No overrepresented sequences found.")

    ## Find the total number of sequences, the average percentage
    ## and number of files it is found in
    df <- dplyr::summarise(
        dplyr::group_by(df),
        Total = sum(Count),
        Percentage = mean(Percentage),
        nFiles = dplyr::n(),
        .groups = "drop"
    )
    df <- dplyr::arrange(df, desc(Percentage))
    df <- dplyr::top_n(df, n, Percentage)
    n <- nrow(df) # In case less than the requested n were present

    ## Create the fasta header & data structure
    hdr <- paste0(
        "> Sequence", seq_len(n),
        " Total:", df$Total,
        " Length:", nchar(df$Sequence),
        " Libraries:", df$nFiles,
        " AvePercentage:", scales::percent(df$Percentage/100)
    )
    fasta <- paste(hdr, df$Sequence, sep = "\n")
    writeLines(fasta, path)
    invisible(fasta)
}
)
