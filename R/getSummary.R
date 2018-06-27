#' Get the summary information from Fastqc Files
#'
#' @description Read the information from the \code{summary.txt} files in each FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList
#'
#' @return A \code{tibble} will be returned when supplying a \code{FastqcFile} object,
#' whilst a list of tibbles will be returned when supplying a \code{FastqcFileList} object
#'
#' @examples
#'
#' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "ngsReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # Return a data_frame/tibble with the raw information
#' getSummary(fdl)
#'
#' @importFrom utils unzip
#'
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "FastqcFile",
          function(object){
            path <- path(object)
            if (isCompressed(path, type = "zip")){
              #Get the internal path within the zip archive
              if (!file.exists(path)) stop("The zip archive can not be found.")
              fl <- file.path( gsub(".zip$", "", fileName(object)), "summary.txt")
              # Check the required file exists
              allFiles <- unzip(path, list = TRUE)$Name
              stopifnot(fl %in% allFiles)
              # Open the connection & read the 12 lines
              uz <- unz(path,fl)
              summaryData <- readLines(uz, 12L)
              close(uz)
              # Form the output
              summaryData <- stringr::str_split_fixed(summaryData, pattern = "\t", n = 3)
              summaryData <- tibble::as_tibble(summaryData)
            }
            else{
              # The existence of this file will have been checked at object instantion
              # Check in case it has been deleted post-instantiation though
              fl <- file.path(path, "summary.txt")
              if (!file.exists(fl)) stop("'summary.txt' could not be found.")
              summaryData <- readr::read_delim(fl, delim = "\t", col_names = FALSE)
            }
            colnames(summaryData) <- c("Status", "Category", "Filename")
            summaryData
          })

#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "FastqcFileList",
          function(object){
            out <- lapply(object, getSummary)
            dplyr::bind_rows(out)
          })

#' @export
#' @rdname getSummary
#' @aliases getSummary
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

#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "FastqcData", function(object){object@Summary})

#' @export
#' @rdname getSummary
#' @aliases getSummary
setMethod("getSummary", "FastqcDataList",
          function(object){
            df <- lapply(object@.Data, getSummary)
            dplyr::bind_rows(df)
          })
