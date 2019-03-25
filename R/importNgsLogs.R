#' @title Import Various NGS-related log files
#'
#' @description Imports NGS-related log files such as those generated from
#' stderr
#'
#' @details Imports one or more log files as output by tools such as:
#' \code{bowtie}
#'
#' @param x \code{character}. Vector of filenames. All log files must be of the
#' same type. Duplicates of the file name will be ignored.
#' @param type \code{character}. The type of file being imported. Can be one of
#' \code{bowtie}
#'
#' @return A \code{tibble}.
#' Column names are broadly similar to the text in supplied files,
#' but have been modified for easier handling under R naming conventions.
#'
#' @examples
#' f <- c("bowtiePE.txt", "bowtieSE.txt")
#' bowtieLogs <- system.file("extdata", f, package = "ngsReports")
#' df <- importNgsLogs(bowtieLogs, type = "bowtie")
#'
#' @importFrom lubridate dminutes
#' @importFrom lubridate dhours
#' @importFrom lubridate dseconds
#'
#' @export
importNgsLogs <- function(x, type) {

    x <- unique(x) # Remove any  duplicates
    stopifnot(file.exists(x)) # Check they all exist

    ## Check for a valid filetype
    possibleTypes <- c("bowtie")
    type <- match.arg(type, possibleTypes)

    ## Load the data
    data <- suppressWarnings(lapply(x, readLines))
    names(data) <- basename(x)

    if (type == "bowtie") {

        ## Perform a check on the data
        validLogs <- vapply(data, .isValidBowtieLog, logical(1))
        if (any(!validLogs)) {
            failed <- names(validLogs)[!validLogs]
            stop(paste("Incorrect file structure for:", failed , collapse = "\n"))
        }

        ## Load the data
        df <- .parseBowtieLogs(data)
    }

    as_tibble(df)

}

#' @title Check for correct structure of supplied Bowtie log files
#' @description Check for correct structure of supplied Bowtie log files after
#' reading in using readLines.
#' @details Checks for all the required fields in the lines provided
#' @param x Character vector as output when readLines to a supplied log file
#' @return logical(1)
#' @keywords internal
.isValidBowtieLog <- function(x){
    nLines <- length(x)
    fields <- c(
        "Time loading forward index",
        "Time loading mirror index:",
        "Seeded quality full-index search:",
        "# reads processed:",
        "# reads with at least one reported alignment:",
        "# reads that failed to align:",
        "Time searching:",
        "Overall time:"
    )
    chk <- vapply(fields, function(f){any(grepl(f, x))}, logical(1))
    all(chk)
}

#' @title Parse data from Bowtie log files
#' @description Parse data from Bowtie log files
#' @details Checks for structure will have been performed
#' @param data List of lines read using readLines on one or more files
#' @return data.frame
#' @keywords internal
.parseBowtieLogs <- function(data){

    ## This will have already been through the validity checks
    df <- lapply(data, function(x){
        x <- gsub("# ", "", x)
        ## Treat the line containing the total reads separately
        total <- grep("Reported .* alignments", x = x, value = TRUE)
        x <- setdiff(x, total)
        ## Split the remaining lines into a nx2 matrix,
        ## then change to titles and remove spaces/dashes
        x <- stringr::str_split_fixed(x, pattern = ": ", 2)
        x[,1] <- stringr::str_to_title(x[,1])
        x[,1] <- gsub("( |-)", "_", x[,1])
        ## Remove those percentages from some of the fields
        x[,2] <- gsub("(.+) \\(.+\\)", "\\1", x[,2])
        ## Return a data.frame
        df <- structure(as.list(x[,2]), names = x[,1])
        df <- as.data.frame(df, stringsAsFactors = FALSE)
    })
    df <- dplyr::bind_rows(df)

    ## Some colnames may have a flag from the original bowtie code
    ## This line will correct that after the title case conversion above
    names(df) <- gsub("__(.)", "_-\\L\\1", names(df), perl = TRUE)

    ## Reformat the columns as integers and durations
    timeCols <- grepl("(Time|Full_Index_Search)", names(df))
    df[!timeCols] <- suppressWarnings(lapply(df[!timeCols], as.integer))
    df[timeCols] <- lapply(df[timeCols], function(x){
        x <- stringr::str_split_fixed(x, ":", 3)
        x <- as.numeric(x)
        x <- matrix(x, ncol = 3)
        dhours(x[,1]) + dminutes(x[,2]) + dseconds(x[,3])
    })
    ## Add the filename, reorder the columns & return a tibble
    df$Filename <- basename(x)
    df <- dplyr::select(
        df,
        "Filename",
        tidyselect::contains("Reads"),
        tidyselect::contains("Time"),
        tidyselect::everything()
    )

    df
}