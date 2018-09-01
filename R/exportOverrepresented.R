#' @title Write fasta of Over-Represented sequences.
#'
#' @description Output overrepresented sequences to disk in fasta format.
#'
#' @details Fasta will contain \code{Filename}, \code{Possible Source}, \code{Percent of total reads}
#'
#' @param x Can be a \code{FastqcData} or \code{FastqcDataList}
#' @param path Path to export the fasta file to. Reverts to a default if not supplied
#' @param n The number of sequences to output
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param noAdapters logical. Remove any sequences identified as possible adapters or primers by FastQC
#'
#' @return Exports to a fasta file, and returns the fasta information invisibly
#'
#' @examples
#' \dontrun{
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Export the top10 Overrepresented Sequences as a single fasta file
#' exportOverrepresented(fdl, path = "top10.fa")
#'
#' }
#'
#' @name exportOverrepresented
#' @rdname exportOverrepresented-methods
#' @export
setGeneric("exportOverrepresented",function(x, path, n, labels, noAdapters = TRUE){standardGeneric("exportOverrepresented")})
#' @aliases exportOverrepresented,FastqcData
#' @rdname exportOverrepresented-methods
#' @export
setMethod("exportOverrepresented", signature = "FastqcData",
          function(x, path,  n = 10, labels, noAdapters = TRUE){
              
              df <- Overrepresented_sequences(x)
              
              if(missing(labels)){
                  labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
              }
              else{
                  if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
              }
              if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
              df$Filename <- labels[df$Filename]
              
              if (missing(path)) path <- paste(unique(df$Filename), "top", n, "overrepresented.fa", sep = "_")
              
              # Remove any putative adapter or primer sequences
              if (noAdapters) df <- df[!grepl("(Primer|Adapter)", df$Possible_Source),]
              if (nrow(df) == 0) stop("No overrepresented sequences found.")
              
              # Now just extract the top n & export
              df <- dplyr::top_n(df, n, Percentage)
              n <- nrow(df) # In case less than the requested n were present
              hdr <- paste0("> ", df$Filename, " Sequence", seq_len(n),
                            " Count:", df$Count,
                            " Length:", nchar(df$Sequence))
              fasta <- paste(hdr, df$Sequence, sep="\n")
              writeLines(fasta, path)
              invisible(fasta)
          }
)
#' @aliases exportOverrepresented,FastqcDataList
#' @rdname exportOverrepresented-methods
#' @export
setMethod("exportOverrepresented", signature = "FastqcDataList",
          function(x, path,  n = 10, labels, noAdapters = TRUE){
              
              df <- Overrepresented_sequences(x)
              
              if(missing(labels)){
                  labels <- structure(
                      gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), 
                      names = unique(df$Filename)
                  )
              }
              else{
                  if (!all(unique(df$Filename) %in% names(labels))) stop(
                      "All file names must be included as names in the vector of labels"
                  )
              }
              if (length(unique(labels)) != length(labels)) stop(
                  "The labels vector cannot contain repeated values"
              )
              df$Filename <- labels[df$Filename]
              
              if (missing(path)) path <- paste(
                  "all_files", "top", n, "overrepresented.fa", sep = "_"
              )
              # Remove any putative adapter or primer sequences
              if (noAdapters) df <- df[!grepl("(Primer|Adapter)", df$Possible_Source),]
              if (nrow(df) == 0) stop("No overrepresented sequences found.")
              
              Count <- Sequence <- c() # Declaration to avoid NOTES during R CMD check
              df <- dplyr::group_by(df, Sequence)
              # Find the total number of sequences, the average percentage and number of files it is found in
              df <- dplyr::summarise(df, Total = sum(Count), 
                                     Percentage = mean(Percentage), nFiles = n())
              df <- dplyr::arrange(df, desc(Percentage))
              df <- dplyr::top_n(df, n, Percentage)
              n <- nrow(df) # In case less than the requested n were present
              
              # Create the fasta header & data structure
              hdr <- paste0("> Sequence", seq_len(n),
                            " Total:", df$Total,
                            " Length:", nchar(df$Sequence),
                            " Libraries:", df$nFiles,
                            " AvePercentage:", scales::percent(df$Percentage/100))
              fasta <- paste(hdr, df$Sequence, sep="\n")
              writeLines(fasta, path)
              invisible(fasta)
          }
)
