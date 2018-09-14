#' @title Perform the checks and return the labels
#'
#' @description Checks for the presence of labels and returns defaults
#'
#' @details Takes a named vector of labels and checks for the correct fields.
#' If no vector is supplied, returns the file names missing the specified pattern,
#' which defaults to removing fastq(.gz), fq(.gz), bam, sam or cram.
#'
#' @param df A data.frame with a column titled "Filename"
#' @param labels Named vector of labels for plotting
#' @param pattern character Regular expression to remove from filenames
#' 
#' @return Named character vector
#'
#' @examples
#' df <- data.frame(Filename = paste0("File1", "File2"), "fastq")
#' ngsReports:::setLabels(df)
#'
#' @keywords internal
setLabels <- function(df, labels, pattern = ".(fastq|fq|bam|sam|cram).*"){
    stopifnot(is.data.frame(df))
    stopifnot("Filename" %in% names(df))
    if (missing(labels)) {
        labels <- structure(gsub(pattern, "", unique(df$Filename)),
                            names = unique(df$Filename))
    }
    if (!all(df$Filename %in% names(labels))) stop(
        "Names of supplied labels must match all filenames."
    )
    if (any(duplicated(labels))) stop(
        "Labels must be unique."
    )
    labels
}
