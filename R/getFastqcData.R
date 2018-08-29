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

            # Define the output to have the same structure as fastqcData, except with
            # an additional slot to account for the changes in the Sequence Duplication Levels output
            out <- vector("list", length = length(modules) + 1)
            names(out) <- c(modules, "Total_Deduplicated_Percentage")

            ## Get the Basic Statistics
            Basic_Statistics <- c()
            if ("Basic_Statistics" %in% modules){
              Basic_Statistics <- getBasicStatistics(fastqcLines)
              stopifnot(is.data.frame(Basic_Statistics))
            }
            out[["Basic_Statistics"]] <- Basic_Statistics

            ## Get the Per Base Sequence Qualities
            Per_base_sequence_quality <- c()
            if ("Per_base_sequence_quality" %in% modules){
              Per_base_sequence_quality <- getPerBaseSeqQuals(fastqcLines)
              stopifnot(is.data.frame(Per_base_sequence_quality))
            }
            out[["Per_base_sequence_quality"]] <- Per_base_sequence_quality

            ## Get the Per Tile Sequence Qualities
            Per_tile_sequence_quality <- c()
            if ("Per_tile_sequence_quality" %in% modules){
              Per_tile_sequence_quality <- getPerTileSeqQuals(fastqcLines)
              stopifnot(is.data.frame(Per_tile_sequence_quality))
            }
            out[["Per_tile_sequence_quality"]] <- Per_tile_sequence_quality

            ## Get the Per Sequence Quality Scores
            Per_sequence_quality_scores <- getPerSeqQualScores(fastqcLines)
            stopifnot(is.data.frame(Per_sequence_quality_scores))
            out[["Per_sequence_quality_scores"]] <- Per_sequence_quality_scores

            ## Get the Per Base Sequence Qualities
            Per_base_sequence_content <- c()
            if ("Per_base_sequence_content" %in% modules){
              Per_base_sequence_content <- getPerBaseSeqContent(fastqcLines)
              stopifnot(is.data.frame(Per_base_sequence_content))
            }
            out[["Per_base_sequence_content"]] <- Per_base_sequence_content

            ## Get the Per Sequence GC Content
            Per_sequence_GC_content <- c()
            if("Per_sequence_GC_content" %in% modules){
              Per_sequence_GC_content <- getPerSeqGcContent(fastqcLines)
              stopifnot(is.data.frame(Per_sequence_GC_content))
            }
            out[["Per_sequence_GC_content"]] <- Per_sequence_GC_content

            ## Get the Per Base Sequence Qualities
            Per_base_N_content <- c()
            if ("Per_base_N_content" %in% modules){
              Per_base_N_content <- getPerBaseNContent(fastqcLines)
              stopifnot(is.data.frame(Per_base_N_content))
            }
            out[["Per_base_N_content"]] <- Per_base_N_content

            ## Get the Sequence Length Distribution
            Sequence_Length_Distribution <- c()
            if ("Sequence_Length_Distribution" %in% modules){
              Sequence_Length_Distribution <- getSeqLengthDist(fastqcLines)
              stopifnot(is.data.frame(Sequence_Length_Distribution))
            }
            out[["Sequence_Length_Distribution"]] <- Sequence_Length_Distribution

            # Get the Sequence Duplication Levels
            Sequence_Duplication_Levels <- c()
            Total_Deduplicated_Percentage <- c()
            if("Sequence_Duplication_Levels" %in% modules){
              Sequence_Duplication_Levels <- getSeqDuplicationLevels(fastqcLines)
              stopifnot(is.data.frame(Sequence_Duplication_Levels[["Sequence_Duplication_Levels"]]))
            }
            out[["Sequence_Duplication_Levels"]] <- Sequence_Duplication_Levels[["Sequence_Duplication_Levels"]]
            out[["Total_Deduplicated_Percentage"]] <- Sequence_Duplication_Levels[["Total_Deduplicated_Percentage"]]

            # Get the Overrepresented Sequences
            Overrepresented_sequences <- c()
            if("Overrepresented_sequences" %in% modules) {
              Overrepresented_sequences <- getOverrepSeq(fastqcLines)
              stopifnot(is.data.frame(Overrepresented_sequences))
            }
            out[["Overrepresented_sequences"]] <- Overrepresented_sequences

            # Get the Adapter Content
            Adapter_Content <- c()
            if("Adapter_Content" %in% modules){
              Adapter_Content <- getAdapterContent(fastqcLines)
              stopifnot(is.data.frame(Adapter_Content))
            }
            out[["Adapter_Content"]] <- Adapter_Content

            # Get the Kmer Content
            Kmer_Content <- c()
            if("Kmer_Content" %in% modules){
              Kmer_Content <- getKmerContent(fastqcLines)
              stopifnot(is.data.frame(Kmer_Content))
            }
            out[["Kmer_Content"]] <- Kmer_Content

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

# Define a series of functions for arranging the data
# after splitting the input from readLines()
getBasicStatistics <- function(fastqcLines){

  x <- fastqcLines[["Basic_Statistics"]][-1]
  vals <- gsub(".+\\t(.+)", "\\1", x)
  names(vals) <- gsub("(.+)\\t.+", "\\1", x)
  names(vals) <- gsub(" ", "_", names(vals))

  # Check for the required values
  reqVals <- c("Filename", "File_type", "Encoding", "Total_Sequences",
               "Sequences_flagged_as_poor_quality", "Sequence_length", "%GC")
  stopifnot(reqVals %in% names(vals))

  df <- tibble::as_tibble(as.list(vals))
  df$Total_Sequences <- as.integer(df$Total_Sequences)
  df$Sequences_flagged_as_poor_quality <- as.integer(df$Sequences_flagged_as_poor_quality)
  df$Shortest_sequence <- as.integer(gsub("(.*)-.*", "\\1", df$Sequence_length))
  df$Longest_sequence <- as.integer(gsub(".*-(.*)", "\\1", df$Sequence_length))
  dplyr::select(df,
                dplyr::one_of("Filename", "Total_Sequences"),
                dplyr::contains("quality"), dplyr::ends_with("sequence"),
                dplyr::one_of("%GC", "File_type", "Encoding"))

}

getPerBaseSeqQuals <- function(fastqcLines){

  x <- fastqcLines[["Per_base_sequence_quality"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Base", "Mean", "Median", "Lower_Quartile", "Upper_Quartile",
               "10th_Percentile", "90th_Percentile")
  stopifnot(reqVals %in% names(df))

  # Change all columns except the position to numeric values
  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
  df

}

getPerTileSeqQuals <- function(fastqcLines){

  x <- fastqcLines[["Per_tile_sequence_quality"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Tile", "Base", "Mean")
  stopifnot(reqVals %in% names(df))

  df[["Mean"]] <- as.numeric(df[["Mean"]])
  df

}

getPerSeqQualScores <- function(fastqcLines){

  x <- fastqcLines[["Per_sequence_quality_scores"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Quality", "Count")
  stopifnot(reqVals %in% names(df))

  df$Quality <- as.integer(df$Quality)
  df$Count <- as.integer(df$Count)
  df
}

getPerBaseSeqContent <- function(fastqcLines){

  x <- fastqcLines[["Per_base_sequence_content"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Base", "G", "A", "T", "C")
  stopifnot(reqVals %in% names(df))

  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
  df

}

getPerSeqGcContent <- function(fastqcLines){

  x <- fastqcLines[["Per_sequence_GC_content"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("GC_Content", "Count")
  stopifnot(reqVals %in% names(df))

  df$GC_Content <- as.integer(df$GC_Content)
  df$Count <- as.integer(df$Count)
  df
}

getPerBaseNContent <- function(fastqcLines){

  x <- fastqcLines[["Per_base_N_content"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Base", "N-Count")
  stopifnot(reqVals %in% names(df))

  df$`N-Count` <- as.integer(df$`N-Count`)
  df

}

getSeqLengthDist <- function(fastqcLines){

  x <- fastqcLines[["Sequence_Length_Distribution"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Length", "Count")
  stopifnot(reqVals %in% names(df))

  df$Lower <- as.integer(gsub("(.*)-.*", "\\1", df$Length))
  df$Upper <- as.integer(gsub(".*-(.*)", "\\1", df$Length))
  df$Count = as.integer(df$Count)

  df[c("Length", "Lower", "Upper", "Count")]

}

getSeqDuplicationLevels <- function(fastqcLines){

  x <- fastqcLines[["Sequence_Duplication_Levels"]][-1]

  # Check for the presence of the Total Deduplicate Percentage value
  hasTotDeDup <- grepl("Total Deduplicated Percentage", x)
  Total_Deduplicated_Percentage <- NA_real_
  if (any(hasTotDeDup)){
    Total_Deduplicated_Percentage <- as.numeric(gsub(".+\\t(.*)", "\\1", x[hasTotDeDup]))
    stopifnot(length(Total_Deduplicated_Percentage) == 1)
  }

  # Remove the Total value entry from the original object
  x <- x[!hasTotDeDup]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Duplication_Level", "Percentage_of_deduplicated", "Percentage_of_total")
  stopifnot(reqVals %in% names(df))

  # Convert to numeric
  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)

  # Return a list with both values
  list(Total_Deduplicated_Percentage = Total_Deduplicated_Percentage,
       Sequence_Duplication_Levels = df)

}

getOverrepSeq <- function(fastqcLines){

  x <- fastqcLines[["Overrepresented_sequences"]][-1]
  if (length(x) <= 1) return(dplyr::data_frame())
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Sequence", "Count", "Percentage", "Possible_Source")
  stopifnot(reqVals %in% names(df))

  df$Count <- as.integer(df$Count)
  df$Percentage <- as.numeric(df$Percentage)
  df
}

getAdapterContent <- function(fastqcLines){

  x <- fastqcLines[["Adapter_Content"]][-1]
  if (length(x) <= 1) return(dplyr::data_frame())
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- "Position"
  stopifnot(reqVals %in% names(df))
  stopifnot(ncol(df) > 1)

  df[!names(df) %in% reqVals] <- lapply(df[!names(df) %in% reqVals], as.numeric)
  df

}

getKmerContent <- function(fastqcLines){

  x <- fastqcLines[["Kmer_Content"]][-1]
  if (length(x) <= 1) return(dplyr::data_frame())
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- tibble::as_tibble(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Sequence", "Count", "PValue", "Obs/Exp_Max", "Max_Obs/Exp_Position")
  stopifnot(reqVals %in% names(df))

  df$Count <- as.integer(df$Count)
  df$PValue <- as.numeric(df$PValue)
  df$`Obs/Exp_Max` <- as.numeric(df$`Obs/Exp_Max`)
  df$`Max_Obs/Exp_Position` <- as.character(df$`Max_Obs/Exp_Position`)

  df

}

# Define a function to split each module by the tab symbol
# By default, this will place the first element as the column names
splitByTab <- function(x, firstRowToNames = TRUE, tab = "\\t"){

  nCol <- stringr::str_count(x[1], pattern = tab) + 1
  if (firstRowToNames){
    # Get the first element as a vector
    nm <- stringr::str_split_fixed(x[1], pattern = tab, n = nCol)
    # Split the remainder
    df <- stringr::str_split_fixed(x[-1], pattern = tab, n = nCol)
    colnames(df) <- nm
  }
  else {
    df <- stringr::str_split_fixed(x, pattern = tab, n = nCol)
  }
  as.data.frame(df)
}
