#' Get all data from FastQC files
#'
#' @description Read the information from the \code{fastqc_data.txt} files in each FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList, or paths to files
#'
#' @return An object of \code{FastqcData} or a \code{FastqcDataList}
#'
#' @include AllGenerics.R
#'
#' @examples
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' @export
#' @rdname getFastqcData-methods
setGeneric("getFastqcData", function(object){standardGeneric("getFastqcData")})

#' @importFrom utils unzip
#' @name getFastqcData
#' @aliases getFastqcData,FastqcFile-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "FastqcFile",
          function(object){

            # Setup the path to the file & use tryCatch to ensure it exists
            # when testing to see if it is compressed
            path <- path(object)
            comp <- tryCatch(isCompressed(path, type = "zip"))

            if (comp){

              # Get the internal path within the zip archive
              fl <- file.path( gsub(".zip$", "", fileName(object)), "fastqc_data.txt")

              # Check the required file exists within the file
              allFiles <- unzip(path, list = TRUE)$Name
              stopifnot(fl %in% allFiles)

              # # Open the connection & read the 12 lines
              uz <- unz(path, fl)
              fastqcLines <- readLines(uz)
              close(uz)

            }
            else{

              # The existence of this file will have been checked at instantiation
              # of the FastqcFile. Check in case it has been deleted post-instantiation
              fl <- file.path(path, "fastqc_data.txt")
              if (!file.exists(fl)) stop("'fastqc_data.txt' could not be found.")
              fastqcLines <- readLines(fl)

            }

            # Modified from Repitools::readFastQC
            # Remove any '#' symbols
            fastqcLines <- gsub("#", "", fastqcLines)
            # Remove any lines which specify '>>END_MODULE'
            fastqcLines <- fastqcLines[!grepl(">>END_MODULE", fastqcLines)]

            # The FastQC version NUMBER
            vers <- strsplit(fastqcLines[[1]], split = "\t")[[1]][2]
            # Remove the version line for easier module identification
            fastqcLines <- fastqcLines[-1]

            # Setup the module names
            modules <- grep("^>>", fastqcLines, value = TRUE)
            # Extract the text after the '>>' & before the tab
            modules <- gsub("^>>(.+)\\t.+", "\\1", modules)
            modules <- gsub(" ", "_", modules) # Add underscores

            # Define the standard modules
            reqModules <- c("Basic_Statistics", "Per_base_sequence_quality",
                            "Per_tile_sequence_quality", "Per_sequence_quality_scores",
                            "Per_base_sequence_content", "Per_sequence_GC_content",
                            "Per_base_N_content", "Sequence_Length_Distribution",
                            "Sequence_Duplication_Levels", "Overrepresented_sequences",
                            "Adapter_Content", "Kmer_Content")
            # Check that at least one of the standard modules is present
            if (!any(modules %in% reqModules)) stop("None of the standard modules were able to be found in the data.")

            # Split the data based on the '>>' pattern, which will indicate the
            # beginning of a new module
            fastqcLines <- split(fastqcLines, cumsum(grepl("^>>", fastqcLines)))
            # Assign the module names
            names(fastqcLines) <- vapply(fastqcLines, function(x){
              # Get the first line before the tab separator
              nm <- gsub(">>(.+)\\t.+", "\\1", x[[1]])
              # Replace the space with underscore
              gsub(" ", "_", nm)
            }, character(1))
            # Remove the name from each module's data
            fastqcLines <- lapply(fastqcLines, function(x){x[-1]})

            # Define the output to have the same structure as fastqcData, except
            # with an additional slot to account for the Sequence Duplication
            # Levels output
            # Initialise an empty list based on the standard modules after checking
            # for the Total_Deduplicated Percentage
            allModules <- unique(c(reqModules, modules, "Total_Deduplicated_Percentage"))
            out <- sapply(allModules, function(x){NULL}, simplify = FALSE)

            ## Get the Modules
            out[["Basic_Statistics"]] <- getBasicStatistics(fastqcLines)
            out[["Per_base_sequence_quality"]] <- getPerBaseSeqQuals(fastqcLines)
            out[["Per_tile_sequence_quality"]] <- getPerTileSeqQuals(fastqcLines)
            out[["Per_sequence_quality_scores"]] <- getPerSeqQualScores(fastqcLines)
            out[["Per_base_sequence_content"]] <- getPerBaseSeqContent(fastqcLines)
            out[["Per_sequence_GC_content"]] <- getPerSeqGcContent(fastqcLines)
            out[["Per_base_N_content"]] <-  getPerBaseNContent(fastqcLines)
            out[["Sequence_Length_Distribution"]] <- getSeqLengthDist(fastqcLines)
            out[["Overrepresented_sequences"]] <- getOverrepSeq(fastqcLines)
            out[["Adapter_Content"]] <- getAdapterContent(fastqcLines)
            out[["Kmer_Content"]] <- getKmerContent(fastqcLines)
            
            # Get the Sequence Duplication Levels
            Sequence_Duplication_Levels <- getSeqDuplicationLevels(fastqcLines)
            dupMods <- names(Sequence_Duplication_Levels)
            out[dupMods] <- Sequence_Duplication_Levels[dupMods]

            #Get the summary
            Summary <- getSummary(object)
            stopifnot(is.data.frame(Summary))

            args <- c(list(Class = "FastqcData",
                           path = path(object),
                           Version = vers,
                           Summary = Summary),
                      out)

            do.call("new", args)

          })

#' @name getFastqcData
#' @aliases getFastqcData,NULL-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "NULL",
          function(object){
            if(is.null(object))stop("No files have been provided, please read in files")
          })

#' @name getFastqcData
#' @aliases getFastqcData,character-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "character",
          function(object){
            if(length(object) ==1) {
              object <- FastqcFile(object)
            }
            else{
              object <- FastqcFileList(object)
            }
            getFastqcData(object)
          })

#' @name getFastqcData
#' @aliases getFastqcData,FastqcFileList-method
#' @rdname getFastqcData-methods
#' @export
setMethod("getFastqcData", "FastqcFileList",
          function(object){
             fqc <- lapply(object@.Data, getFastqcData)
             new("FastqcDataList", fqc)
          })

# Define a series of quick functions for arranging the data
# after splitting the input from readLines()
getBasicStatistics <- function(fastqcLines){

  # Return NULL if the module is missing
  if (!"Basic_Statistics" %in% names(fastqcLines)) return(NULL)

  x <- splitByTab(fastqcLines[["Basic_Statistics"]])
  stopifnot(setequal(names(x), c("Measure", "Value")))
  vals <- x[["Value"]]
  names(vals) <- gsub(" ", "_", x[["Measure"]])

  # Check for the required values
  reqVals <- c("Filename", "File_type", "Encoding", "Total_Sequences",
               "Sequences_flagged_as_poor_quality", "Sequence_length", "%GC")
  stopifnot(reqVals %in% names(vals))

  # Setup the data_frame with the correct value types
  df <- tibble::as_tibble(as.list(vals))
  df$Shortest_sequence <- gsub("(.*)-.*", "\\1", df$Sequence_length)
  df$Longest_sequence <- gsub(".*-(.*)", "\\1", df$Sequence_length)
  intVals <- c("Shortest_sequence", "Longest_sequence", "Total_Sequences",
               "Sequences_flagged_as_poor_quality")
  df[intVals] <-lapply(df[intVals], as.integer)
  dplyr::select(df,
                "Filename", "Total_Sequences",
                dplyr::contains("quality"), dplyr::ends_with("sequence"),
                dplyr::one_of("%GC", "File_type", "Encoding"))

}

getPerBaseSeqQuals <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Per_base_sequence_quality"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])
  names(df) <- gsub(" ", "_", names(df))

  # Check for the required values
  reqVals <- c("Base", "Mean", "Median", "Lower_Quartile", "Upper_Quartile",
               "10th_Percentile", "90th_Percentile")
  stopifnot(reqVals %in% names(df))

  # Change to numeric where appropriate
  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
  tibble::as_tibble(df)

}

getPerTileSeqQuals <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Per_tile_sequence_quality"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])

  # Check for the required values
  reqVals <- c("Tile", "Base", "Mean")
  stopifnot(reqVals %in% names(df))

  df[["Mean"]] <- as.numeric(df[["Mean"]])
  tibble::as_tibble(df)

}

getPerSeqQualScores <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Per_sequence_quality_scores"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])

  # Check for the required values
  reqVals <- c("Quality", "Count")
  stopifnot(reqVals %in% names(df))

  df <- lapply(df, as.integer)
  tibble::as_tibble(df)
}

getPerBaseSeqContent <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Per_base_sequence_content"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])

  # Check for the required values
  reqVals <- c("Base", "G", "A", "T", "C")
  stopifnot(reqVals %in% names(df))

  # Set all as numeric except the 'Base'column
  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
  tibble::as_tibble(df)

}

getPerSeqGcContent <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Per_sequence_GC_content"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])
  names(df) <- gsub(" ", "_", names(df))

  # Check for the required values
  reqVals <- c("GC_Content", "Count")
  stopifnot(reqVals %in% names(df))

  # Convert to integers
  df[["GC_Content"]] <- as.integer(df[["GC_Content"]])
  df[["Count"]] <- as.numeric(df[["Count"]])
  tibble::as_tibble(df)
}

getPerBaseNContent <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Per_base_N_content"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])

  # Check for the required values
  reqVals <- c("Base", "N-Count")
  stopifnot(reqVals %in% names(df))

  df[["N-Count"]] <- as.integer(df[["N-Count"]])
  tibble::as_tibble(df)

}

getSeqLengthDist <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Sequence_Length_Distribution"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)

  df <- splitByTab(fastqcLines[[mod]])

  # Check for the required values
  reqVals <- c("Length", "Count")
  stopifnot(reqVals %in% names(df))

  # If there is a range, separate into Lower & Upper values
  df$Lower <- as.integer(gsub("(.*)-.*", "\\1", df$Length))
  df$Upper <- as.integer(gsub(".*-(.*)", "\\1", df$Length))
  df$Count <- as.integer(df$Count)
  df <- tibble::as_tibble(df)

  df[c("Length", "Lower", "Upper", "Count")]

}

getSeqDuplicationLevels <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Sequence_Duplication_Levels"
  if (!mod %in% names(fastqcLines)) {
    return(list(Total_Deduplicated_Percentage = NULL,
                Sequence_Duplication_Levels = NULL))
  }
  if (length(fastqcLines[[mod]]) == 0) {
    return(list(Total_Deduplicated_Percentage = NULL,
                Sequence_Duplication_Levels = NULL))
  }
  x <- fastqcLines[[mod]]

  # Check for the presence of the Total Deduplicate Percentage value
  hasTotDeDup <- grepl("Total Deduplicated Percentage", x)
  Total_Deduplicated_Percentage <- NA_real_
  if (any(hasTotDeDup)){
    Total_Deduplicated_Percentage <- as.numeric(gsub(".+\\t(.*)", "\\1", x[hasTotDeDup]))
  }
  if (length(Total_Deduplicated_Percentage) > 1) stop("Too many elements matched Total_Deduplicated_Percentage")

  # Remove the Total value entry from the original object
  df <- splitByTab(x[!hasTotDeDup])
  names(df) <- gsub(" ", "_", names(df))

  # Check for the required values
  reqVals <- c("Duplication_Level", "Percentage_of_deduplicated", "Percentage_of_total")
  stopifnot(reqVals %in% names(df))

  # Convert percentages to numeric
  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)

  # Return a list with both values
  list(Total_Deduplicated_Percentage = Total_Deduplicated_Percentage,
       Sequence_Duplication_Levels = tibble::as_tibble(df))

}

getOverrepSeq <- function(fastqcLines){

  # Return NULL if the module is missing
  mod <- "Overrepresented_sequences"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)
  
  df <- splitByTab(fastqcLines[[mod]])
  names(df) <- gsub(" ", "_", names(df))

  # Check for the required values
  reqVals <- c("Sequence", "Count", "Percentage", "Possible_Source")
  stopifnot(reqVals %in% names(df))

  df$Count <- as.integer(df$Count)
  df$Percentage <- as.numeric(df$Percentage)
  tibble::as_tibble(df)

}

getAdapterContent <- function(fastqcLines){
  
  # Return NULL if the module is missing
  mod <- "Adapter_Content"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)

  df <- splitByTab(fastqcLines[[mod]])
  names(df) <- gsub(" ", "_", names(df))

  # Check for the required values
  reqVals <- "Position"
  stopifnot(reqVals %in% names(df))
  stopifnot(ncol(df) > 1)

  numCols <- !names(df) %in% reqVals
  df[numCols] <- lapply(df[numCols], as.numeric)
  tibble::as_tibble(df)

}

getKmerContent <- function(fastqcLines){
  
  # Return NULL if the module is missing
  mod <- "Kmer_Content"
  if (!mod %in% names(fastqcLines)) return(NULL)
  if (length(fastqcLines[[mod]]) == 0) return(NULL)

  df <- splitByTab(fastqcLines[[mod]])
  names(df) <- gsub(" ", "_", names(df))

  # Check for the required values
  reqVals <- c("Sequence", "Count", "PValue", "Obs/Exp_Max", "Max_Obs/Exp_Position")
  stopifnot(reqVals %in% names(df))

  df$Count <- as.integer(df$Count)
  df$PValue <- as.numeric(df$PValue)
  df$`Obs/Exp_Max` <- as.numeric(df$`Obs/Exp_Max`)

  tibble::as_tibble(df)

}

