#' Get the complete information from FastQC files
#'
#' @description Read the information from the \code{fastqc_data.txt} files in each FastqcFile
#'
#' @param object Can be a FastqcFile or FastqcFileList.
#'
#' @return An object of \code{FastqcData}
#'
#' @importFrom dplyr data_frame
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_split
#'
#' @include AllGenerics.R
#'
#' @export
#' @rdname getFastqcData
setMethod("getFastqcData", "FastqcFile",
          function(object){
            if (isCompressed(object)){
              # Get the internal path within the zip archive
              if (!file.exists(path(object))) stop("The zip archive can not be found.")
              fl <- file.path( gsub(".zip$", "", fileNames(object)), "fastqc_data.txt")
              # Check the required file exists
              allFiles <- unzip(path(object), list = TRUE)$Name
              stopifnot(fl %in% allFiles)
              # # Open the connection & read the 12 lines
              uz <- unz(path(object),fl)
              fastqcData <- readLines(uz)
              close(uz)
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

            # Check all required modules are present in the data
            reqModules <- c("Basic_Statistics", "Per_base_sequence_quality",
                            "Per_tile_sequence_quality", "Per_sequence_quality_scores",
                            "Per_base_sequence_content", "Per_sequence_GC_content",
                            "Per_base_N_content", "Sequence_Length_Distribution",
                            "Sequence_Duplication_Levels", "Overrepresented_sequences",
                            "Adapter_Content", "Kmer_Content")
            stopifnot(reqModules %in% modules)

            # Split the data
            fastqcData <- split(fastqcData, cumsum(grepl("^>>", fastqcData)))[-1]
            names(fastqcData) <- modules

            # Define the output to have the same structure as fastqcData, except with
            # an additional slot to account for the changes in the Sequence Duplication Levels output
            out <- vector("list", length = length(modules) + 1)
            names(out) <- c(modules, "Total_Deduplicated_Percentage")

            ## Get the Basic Statistics
            Basic_Statistics <- getBasicStatistics(fastqcData)
            stopifnot(is.data.frame(Basic_Statistics))
            out[["Basic_Statistics"]] <- Basic_Statistics

            ## Get the Per Base Sequence Qualities
            Per_base_sequence_quality <- getPerBaseSeqQuals(fastqcData)
            stopifnot(is.data.frame(Per_base_sequence_quality))
            out[["Per_base_sequence_quality"]] <- Per_base_sequence_quality

            ## Get the Per Tile Sequence Qualities
            Per_tile_sequence_quality <- getPerTileSeqQuals(fastqcData)
            stopifnot(is.data.frame(Per_tile_sequence_quality))
            out[["Per_tile_sequence_quality"]] <- Per_tile_sequence_quality

            ## Get the Per Sequence Quality Scores
            Per_sequence_quality_scores <- getPerSeqQualScores(fastqcData)
            stopifnot(is.data.frame(Per_sequence_quality_scores))
            out[["Per_sequence_quality_scores"]] <- Per_sequence_quality_scores

            ## Get the Per Base Sequence Qualities
            Per_base_sequence_content <- getPerBaseSeqContent(fastqcData)
            stopifnot(is.data.frame(Per_base_sequence_content))
            out[["Per_base_sequence_content"]] <- Per_base_sequence_content

            ## Get the Per Sequence GC Content
            Per_sequence_GC_content <- getPerSeqGcContent(fastqcData)
            stopifnot(is.data.frame(Per_sequence_GC_content))
            out[["Per_sequence_GC_content"]] <- Per_sequence_GC_content

            ## Get the Per Base Sequence Qualities
            Per_base_N_content <- getPerBaseNContent(fastqcData)
            stopifnot(is.data.frame(Per_base_N_content))
            out[["Per_base_N_content"]] <- Per_base_N_content

            ## Get the Sequence Length Distribution
            Sequence_Length_Distribution <- getSeqLengthDist(fastqcData)
            stopifnot(is.data.frame(Sequence_Length_Distribution))
            out[["Sequence_Length_Distribution"]] <- Sequence_Length_Distribution

            # Get the Sequence Duplication Levels
            Sequence_Duplication_Levels <- getSeqDuplicationLevels(fastqcData)
            stopifnot(is.data.frame(Sequence_Duplication_Levels[["Sequence_Duplication_Levels"]]))
            out[["Sequence_Duplication_Levels"]] <- Sequence_Duplication_Levels[["Sequence_Duplication_Levels"]]
            out[["Total_Deduplicated_Percentage"]] <- Sequence_Duplication_Levels[["Total_Deduplicated_Percentage"]]

            # Get the Overrepresented Sequences
            Overrepresented_sequences <- getOverrepSeq(fastqcData)
            stopifnot(is.data.frame(Overrepresented_sequences))
            out[["Overrepresented_sequences"]] <- Overrepresented_sequences

            # Get the Adapter Content
            Adapter_Content <- getAdapterContent(fastqcData)
            stopifnot(is.data.frame(Adapter_Content))
            out[["Adapter_Content"]] <- Adapter_Content

            # Get the Kmer Content
            Kmer_Content <- getKmerContent(fastqcData)
            stopifnot(is.data.frame(Kmer_Content))
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


#' @export
# #' @rdname getFastqcData
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

#' @export
# #' @rdname getFastqcData
setMethod("getFastqcData", "FastqcFileList",
          function(object){
             fqc <- lapply(object@.Data, getFastqcData)
             new("FastqcDataList", fqc)
          })

# Define a series of functions for arranging the data
# after splitting the input from readLines()
getBasicStatistics <- function(fastqcData){

  x <- fastqcData[["Basic_Statistics"]][-(1:2)]
  vals <- gsub(".+\\t(.+)", "\\1", x)
  names(vals) <- gsub("(.+)\\t.+", "\\1", x)
  names(vals) <- gsub(" ", "_", names(vals))

  # Check for the required values
  reqVals <- c("Filename", "File_type", "Encoding", "Total_Sequences",
               "Sequences_flagged_as_poor_quality", "Sequence_length", "%GC")
  stopifnot(reqVals %in% names(vals))

  df <- dplyr::as_data_frame(as.list(vals))
  df <- dplyr::mutate(df,
                      Total_Sequences = as.integer(Total_Sequences),
                      Sequences_flagged_as_poor_quality = as.integer(Sequences_flagged_as_poor_quality),
                      Shortest_sequence = as.integer(gsub("(.*)-.*", "\\1", Sequence_length)),
                      Longest_sequence = as.integer(gsub(".*-(.*)", "\\1", Sequence_length)))
  dplyr::select(df,
                Filename, Total_Sequences,
                dplyr::contains("quality"), dplyr::ends_with("sequence"),
                `%GC`, File_type, Encoding)

}

getPerBaseSeqQuals <- function(fastqcData){

  x <- fastqcData[["Per_base_sequence_quality"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Base", "Mean", "Median", "Lower_Quartile", "Upper_Quartile",
               "10th_Percentile", "90th_Percentile")
  stopifnot(reqVals %in% names(df))

  # Change all columns except the positionto numeric values
  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
  df

}

getPerTileSeqQuals <- function(fastqcData){

  x <- fastqcData[["Per_tile_sequence_quality"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Tile", "Base", "Mean")
  stopifnot(reqVals %in% names(df))

  df[["Mean"]] <- as.numeric(df[["Mean"]])
  df

}

getPerSeqQualScores <- function(fastqcData){

  x <- fastqcData[["Per_sequence_quality_scores"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Quality", "Count")
  stopifnot(reqVals %in% names(df))

  dplyr::mutate(df, Quality = as.integer(Quality), Count = as.integer(Count))

}

getPerBaseSeqContent <- function(fastqcData){

  x <- fastqcData[["Per_base_sequence_content"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Base", "G", "A", "T", "C")
  stopifnot(reqVals %in% names(df))

  df[reqVals[-1]] <- lapply(df[reqVals[-1]], as.numeric)
  df

}

getPerSeqGcContent <- function(fastqcData){

  x <- fastqcData[["Per_sequence_GC_content"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("GC_Content", "Count")
  stopifnot(reqVals %in% names(df))

  dplyr::mutate(df, GC_Content = as.integer(GC_Content), Count = as.numeric(Count))

}

getPerBaseNContent <- function(fastqcData){

  x <- fastqcData[["Per_base_N_content"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Base", "N-Count")
  stopifnot(reqVals %in% names(df))

  dplyr::mutate(df, `N-Count` = as.numeric(`N-Count`))

}

getSeqLengthDist <- function(fastqcData){

  x <- fastqcData[["Sequence_Length_Distribution"]][-1]
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Length", "Count")
  stopifnot(reqVals %in% names(df))

  df <- dplyr::mutate(df,
                      Lower = as.integer(gsub("(.*)-.*", "\\1",Length)),
                      Upper = as.integer(gsub(".*-(.*)", "\\1",Length)),
                      Count = as.integer(Count))
  dplyr::select(df, Length, Lower, Upper, Count)
}

getSeqDuplicationLevels <- function(fastqcData){

  x <- fastqcData[["Sequence_Duplication_Levels"]][-1]

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
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
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

getOverrepSeq <- function(fastqcData){

  x <- fastqcData[["Overrepresented_sequences"]][-1]
  if (length(x) <= 1) return(dplyr::data_frame())
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Sequence", "Count", "Percentage", "Possible_Source")
  stopifnot(reqVals %in% names(df))

  dplyr::mutate(df,
                Count = as.integer(Count),
                Percentage = as.numeric(Percentage))
}

getAdapterContent <- function(fastqcData){

  x <- fastqcData[["Adapter_Content"]][-1]
  if (length(x) <= 1) return(dplyr::data_frame())
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- "Position"
  stopifnot(reqVals %in% names(df))
  stopifnot(ncol(df) > 1)

  df[!names(df) %in% reqVals] <- lapply(df[!names(df) %in% reqVals], as.numeric)
  df

}

getKmerContent <- function(fastqcData){

  x <- fastqcData[["Kmer_Content"]][-1]
  if (length(x) <= 1) return(dplyr::data_frame())
  mat <- stringr::str_split(x, pattern = "\t", simplify = TRUE)
  nc <- ncol(mat)
  df <- dplyr::as_data_frame(matrix(mat[-1,], ncol= nc))
  names(df) <- gsub(" ", "_", mat[1,])

  # Check for the required values
  reqVals <- c("Sequence", "Count", "PValue", "Obs/Exp_Max", "Max_Obs/Exp_Position")
  stopifnot(reqVals %in% names(df))

  dplyr::mutate(df,
                Count = as.integer(Count),
                PValue = as.numeric(PValue),
                `Obs/Exp_Max` = as.numeric(`Obs/Exp_Max`),
                `Max_Obs/Exp_Position` = as.integer(`Max_Obs/Exp_Position`))

}


