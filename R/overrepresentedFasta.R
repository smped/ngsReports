#' @title Write fasta of Over-Represented sequences.
#'
#' @description Output overrepresented sequences to disk in fasta format.
#'
#' @details Fasta will contain \code{Filename}, \code{Possible Source}, \code{Percent of total reads}
#'
#' @param x Can be a \code{FastqcData} or \code{FastqcDataList}
#' @param outDir Directory to output fasta files to
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param n The number of sequences to plot from an individual file
#' @param ... set generic args
#'
#' @return A standard fasta written to \code{outDir}
#'
#'
#' @importFrom dplyr top_n
#'
#' @name overrepresentedFasta
#' @rdname overrepresentedFasta-methods
#' @export
setGeneric("overrepresentedFasta",function(x, outDir = "./", ...){standardGeneric("overrepresentedFasta")})
#' @aliases overrepresentedFasta,FastqcData
#' @rdname overrepresentedFasta-methods
#' @export
setMethod("overrepresentedFasta", signature = "FastqcData",
          function(x, outDir = "./",  n = 10, labels, ...){

            if(base::missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            df <- Overrepresented_sequences(x)
            df$Filename <- labels[df$Filename]
            df <- top_n(df, n, Percentage)
            df["fastaNames"] <- paste0(">", df$Filename)
            df["fastaNames"] <- paste(df$fastaNames, df$Possible_Source, round(df$Percentage, 4), sep = ", ")
            fasta <- df[c("fastaNames", "Sequence")]

            write.table(fasta, paste(outDir, df$Filename[1], "_OverRepresentedSequences.fa"), col.names = FALSE,
                        row.names = FALSE, quote = FALSE, sep = "\n")
            print(paste("output fasta for", df$Filename[1], "to", outDir))
          }
)
#' @aliases overrepresentedFasta,FastqcDataList
#' @rdname overrepresentedFasta-methods
#' @export
setMethod("overrepresentedFasta", signature = "FastqcDataList",
          function(x, outDir = "./", n = 10, labels, ...){

            if(base::missing(labels)){
              labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
            }
            else{
              if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
            }
            if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

            df <- Overrepresented_sequences(fdl)
            df$Filename <- labels[df$Filename]

            out <- lapply(split(df, df$Filename), function(x){
              x <- top_n(x, n, Percentage)
              x["fastaNames"] <- paste0(">", x$Filename)
              x["fastaNames"] <- paste(x$fastaNames, x$Possible_Source, round(x$Percentage, 4), sep = ", ")
              fasta <- x[c("fastaNames", "Sequence")]

              write.table(fasta, paste(outDir, x$Filename[1], "_OverRepresentedSequences.fa"), col.names = FALSE,
                          row.names = FALSE, quote = FALSE, sep = "\n")
              print(paste("output fasta for", x$Filename[1], "to", outDir))
          }
          )
          })
