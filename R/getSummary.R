#' Get the summary information from Fastqc Files
#'
#' @description Read the information from the \code{summary.txt} files in each FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList
#'
#' @return A \code{tibble} will be returned when supplying a \code{FastqcFile} object,
#' whilst a list of tibbles will be returned when supplying a \code{FastqcFileList} object
#'
#' @importFrom readr read_delim
#' @importFrom stringr str_split_fixed
#' @importFrom tibble as_tibble
#'
#' @include AllGenerics.R
#'
#' @export
#' @rdname getSummary
#' @aliases getSummary,FastqcFile-method
setMethod("getSummary", "FastqcFile",
          function(object){
            if (isCompressed(object)){
              #Get the internal path within the zip archive
              if (!file.exists(path(object))) stop("The zip archive can not be found.")
              fl <- file.path( gsub(".zip$", "", fileNames(object)), "summary.txt")
              # Check the required file exists
              allFiles <- unzip(path(object), list = TRUE)$Name
              stopifnot(fl %in% allFiles)
              # Open the connection & read the 12 lines
              uz <- unz(path(object),fl)
              summaryData <- readLines(uz, 12L)
              close(uz)
              # Form the output
              summaryData <- stringr::str_split_fixed(summaryData, pattern = "\t", n = 3)
              summaryData <- tibble::as_tibble(summaryData)
            }
            else{
              # The existence of this file will have been checked at object instantion
              # Check in case it has been deleted post-instantiation though
              fl <- file.path(path(object), "summary.txt")
              if (!file.exists(fl)) stop("'summary.txt' could not be found.")
              summaryData <- readr::read_delim(fl, delim = "\t", col_names = FALSE)
            }
            colnames(summaryData) <- c("Status", "Category", "Filename")
            summaryData
          })

#' @export
#' @rdname getSummary
#' @aliases getSummary,FastqcFileList-method
setMethod("getSummary", "FastqcFileList",
          function(object){
            out <- lapply(object, getSummary)
            dplyr::bind_rows(out)
          })

#' @export
setMethod("getSummary", "character",
          function(object){
            if(length(object) ==1) {
              object <- FastqcFile(object)
            }
            else{
              object <- FastqcFileList(object)
            }
            getSummary(object)
          })
