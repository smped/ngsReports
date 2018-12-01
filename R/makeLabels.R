#' @title Perform the checks and return the labels
#'
#' @description Checks for the presence of labels and returns defaults
#'
#' @details Takes a named vector of labels and checks for the correct fields.
#' If no vector is supplied, returns the file names missing the specified
#' pattern, which defaults to removing the suffixes fastq(.gz), fq(.gz),
#' bam, sam or cram.
#'
#' @param df A data.frame with a column titled "Filename"
#' @param labels Named vector of labels for plotting
#' @param pattern character Regular expression to remove from filenames
#' @param col character Column to use for generating labels
#' @param ... Not used
#'
#' @return Named character vector
#'
#' @examples
#' f <- paste0(c("File1", "File2"), ".fastq")
#' df <- data.frame(Filename = f, stringsAsFactors = FALSE)
#' ngsReports:::makeLabels(df)
#'
#' @keywords internal
makeLabels <- function(df, labels, pattern = ".(fastq|fq|bam|sam|cram).*",
                      col ="Filename", ...){
    stopifnot(is.data.frame(df))
    col <- match.arg(col, colnames(df))
    ## If no labels are provided, just remove the file suffix as determined by
    ## the supplied pattern
    if (missing(labels)) {
        labels <- structure(gsub(pattern, "", unique(df[[col]])),
                            names = unique(df[[col]]))
    }
    if (!all(df[[col]] %in% names(labels))) stop(
        "Names of supplied labels must match all filenames."
    )
    if (any(duplicated(labels))) stop(
        "Labels must be unique."
    )
    labels
}
