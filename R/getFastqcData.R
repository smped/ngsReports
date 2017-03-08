#' Get the complete information from FastQC files
#'
#' @description Read the information from the \code{fastqc_data.txt} files in each FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList
#'
#' @return A \code{tibble} will be returned when supplying a \code{FastqcFile} object,
#' whilst a list of tibbles will be returned when supplying a \code{FastqcFileList} object
#'
#' @importFrom dplyr data_frame
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_split
#'
#' @export
#' @rdname getFastqcData
#' @aliases getFastqcData,FastqcFile-method
setMethod("getFastqcData", "FastqcFile",
          function(object){
            if (isCompressed(object)){
              # #Get the internal path within the zip archive
              # if (!file.exists(path(object))) stop("The zip archive can not be found.")
              # fl <- file.path( gsub(".zip$", "", names(object)), "fastqc_data.txt")
              # # Open the connection & read the 12 lines
              # uz <- unz(path(object),fl)
              # fastqcData <- readLines(uz)
              # close(uz)
              # # Form the output
              # fastqcData <- stringr::str_split_fixed(fastqcData, pattern = "\t", n = 3)
              # fastqcData <- tibble::as_tibble(fastqcData)
            }
            else{
              # The existence of this file will have been checked at object instantion
              # Check in case it has been deleted post-instantiation though
              fl <- file.path(path(object), "fastqc_data.txt")
              if (!file.exists(fl)) stop("'fastqc_data.txt' could not be found.")
              fastqcData <- readLines(fl)
            }
            # Modified from Repitools::readFastQC
            fastqcData <- gsub("#", "", fastqcData)
            fastqcData <- fastqcData[!grepl(">>END_MODULE", fastqcData)]

            # The FastQC version NUMBER
            vers <- gsub("\\t", "", fastqcData[1])

            # Setup the module names
            modules <- grep("^>>", fastqcData, value = TRUE)
            modules <- gsub("^>>", "", modules) # Remove the leading '>>'
            modules <- gsub("(.+)\\t.+", "\\1", modules) # Grab the information before the \t
            modules <- gsub(" ", "_", modules) # Add underscores

            # Check the modules are all present as required

            # Split the data
            fastqcData <- split(fastqcData, cumsum(grepl("^>>", fastqcData)))[-1]
            names(fastqcData) <- modules
            # Define the output to have the same structure as fastqcData, except with
            # an additional slot to account for the changes in the Sequence Duplication Levels output
            out <- vector("list", length = length(modules) + 1)
            names(out) <- c(modules, "Total_Deduplicated_Percentage")

            # Basic Statistics & Sequence_Length_Distribution need to be handled separately
            # This structure is slightly changed from that as provided by FastQC output
            # Sequence length will be separated into `Shortest Sequence` and `Longest Sequence`
            if ("Basic_Statistics" %in% modules){
              nm <- unlist(stringr::str_split(fastqcData$Basic_Statistics[2], pattern = "\t"))
              mat <- stringr::str_split_fixed(fastqcData$Basic_Statistics[-(1:2)], pattern = "\t", n = length(nm))
              v <- mat[,2]
              names(v) <- mat[,1]
              out$Basic_Statistics <- dplyr::as_data_frame(as.list(v))
              out$Basic_Statistics <- dplyr::mutate(out$Basic_Statistics,
                                                    `Total Sequences` = as.integer(`Total Sequences`),
                                                    `Sequences flagged as poor quality` = as.integer(`Sequences flagged as poor quality`),
                                                    `Shortest Sequence` = as.integer(gsub("(.*)-.*", "\\1", `Sequence length`)),
                                                    `Longest Sequence` = as.integer(gsub(".*-(.*)", "\\1", `Sequence length`)))
              out$Basic_Statistics <- dplyr::select(out$Basic_Statistics, -`Sequence length`)
            }
            if ("Sequence_Duplication_Levels" %in% modules){
              tdp <- grepl("Total Deduplicated Percentage", fastqcData$Sequence_Duplication_Levels)
              if (any(tdp)) {
                out$Total_Deduplicated_Percentage <- as.numeric(gsub("Total Deduplicated Percentage\\t(.+)", "\\1",
                                                                    fastqcData$Sequence_Duplication_Levels[tdp] ))
              }
              nm <- unlist(stringr::str_split(fastqcData$Sequence_Duplication_Levels[!tdp][2], pattern = "\t"))
              out$Sequence_Duplication_Levels <- stringr::str_split_fixed(fastqcData$Sequence_Duplication_Levels[!tdp][-(1:2)],
                                                                          pattern = "\t", n = length(nm))
              colnames(out$Sequence_Duplication_Levels) <- nm
              out$Sequence_Duplication_Levels <- dplyr::as_data_frame(out$Sequence_Duplication_Levels)
              out$Sequence_Duplication_Levels$`Percentage of deduplicated` <- as.numeric(out$Sequence_Duplication_Levels$`Percentage of deduplicated`)
              out$Sequence_Duplication_Levels$`Percentage of total` <- as.numeric(out$Sequence_Duplication_Levels$`Percentage of total`)
            }

            # RUN THE REST OF THE MODULES
            # sapply(setdiff(modules, "Sequence_Duplication_Levels"), function(x){
            #   if (length(fastqcData[[x]]) == 1) return(data.frame())
            #   nm <- unlist(stringr::str_split(fastqcData[[x]][2], pattern = "\t"))
            #   nm
            #   # df <- strsplit(fastqcData[[x]][-1], split = "\t")
            # }, simplify = FALSE)
          })
