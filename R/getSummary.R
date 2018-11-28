#' Get the summary information from Fastqc Files
#'
#' @description Read the information from the \code{summary.txt} files in each 
#' FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList
#'
#' @return A \code{tibble} will be returned when supplying a \code{FastqcFile} 
#' object, whilst a list of tibbles will be returned when supplying a 
#' \code{FastqcFileList} object
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
#' getSummary(fdl)
#'
#' @importFrom utils unzip
#' @import tibble
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
              modules <- c("Basic Statistics", 
                           "Per base sequence quality",
                           "Per tile sequence quality", 
                           "Per sequence quality scores",
                           "Per base sequence content", 
                           "Per sequence GC content",
                           "Per base N content", 
                           "Sequence Length Distribution",
                           "Sequence Duplication Levels", 
                           "Overrepresented sequences",
                           "Adapter Content", 
                           "Kmer Content")
              path <- path(object)
              if (isCompressed(path, type = "zip")){
                  #Get the internal path within the zip archive
                  if (!file.exists(path)) stop(
                      "The zip archive can not be found."
                      )
                  fl <- file.path( gsub(".zip$", "", fileName(object)),
                                   "summary.txt")
                  # Check the required file exists
                  allFiles <- unzip(path, list = TRUE)$Name
                  if(!fl %in% allFiles) stop(
                      "summary.txt is missing from the zip archive"
                  )
                  # Open the connection & read all lines
                  uz <- unz(path,fl)
                  summaryData <- readLines(uz)
                  close(uz)
                  
                  # Form the output
                  summaryData <- stringr::str_split_fixed(string = summaryData, 
                                                          pattern = "\t", 
                                                          n = 3)
                  summaryData <- as_tibble(summaryData)
              }
              else{
                  # The existence of this file will have been checked at object
                  # instantion
                  # Check in case it has been deleted post-instantiation though
                  fl <- file.path(path, "summary.txt")
                  if (!file.exists(fl)) stop("'summary.txt' could not be found.")
                  summaryData <- readr::read_delim(file = fl, 
                                                   delim = "\t", 
                                                   col_names = FALSE)
              }
              colnames(summaryData) <- c("Status", "Category", "Filename")
              if (!any(modules %in% summaryData$Category)) stop(
                  "summary.txt contains incomplete data. No expected modules were found."
              )
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
                  tryCatch(object <- FastqcFile(object))
              }
              else{
                  tryCatch(object <- FastqcFileList(object))
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
