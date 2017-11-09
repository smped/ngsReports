#' @title Plot a summary of OVer-represented Sequences
#'
#' @description Plot a summary of Over-represented Sequences for a set of FASTQC reports
#'
#' @details Percentages are obtained by simply summing those within a report.
#' Any possible double counting by FastQC is ignored for the purposes of a simple approximation.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param sequenceSource A regular expression used to group the possible overrepresented sequence source into those
#' matching the pattern, and "Other".
#' Those matching will be classified using the pattern itself, without any supplied braces.
#' Defaults to \code{sequenceSource = "(Primer|Adapter)"}
#' @param col1 The colour to use for the 'Other' percentage
#' @param col2 The colour to use for the percentage matching \code{sequenceSource}
#'
#' @return A standard ggplot2 object
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
#' # Another example which isn't ideal
#' plotOverrepresentedSummary(fdl)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 coord_flip
#' @importFrom grDevices rgb
#'
#' @export
plotOverrepresentedSummary <- function(x, usePlotly = FALSE, labels, sequenceSource = "(Primer|Adapter)",
                                       col1 = rgb(0.2, 0.2, 0.8), col2 = rgb(0.9, 0.2, 0.2)){

  df <- tryCatch(Overrepresented_sequences(x))

  if (nrow(df) == 0) stop("No overrepresented sequences were detected by FastQC")

  # Drop the suffix, or check the alternate labels
  if(base::missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", unique(df$Filename)), names = unique(df$Filename))
  }
  else{
    if (!all(unique(df$Filename) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")


  Percentage <- c() # Here to avoid a NOTE in R CMD Check
  df$Type <- c("Other", sequenceSource)[grepl(sequenceSource, df$Possible_Source) + 1]
  df$Type <- gsub("(\\(|\\))", "", df$Type) # Remove brackets
  df <- dplyr::group_by(df, Filename, Type)
  df <- dplyr::summarise(df, Percentage = sum(Percentage))
  df$Filename <- labels[df$Filename]
  df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

  # Set the axis limits. Just scale the upper limit by 1.05
  ymax <- max(dplyr::summarise(dplyr::group_by(df, Filename), Total = sum(Percentage))$Total)*1.05

  overPlot <- ggplot(df, aes_string(x = "Filename", y = "Percentage", fill = "Type")) +
    geom_bar(stat = "identity") +
    ylab("Overrepresented Sequences (% of Total)") +
    scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    scale_fill_manual(values = c(col1, col2)) +
    theme_bw() +
    coord_flip()

  if (usePlotly){
    message("usePlotly is not yet implemented")

    ###################################################################################
    # This looks terrible. I don't think it can handle the coord flip or stacked bars #
    ###################################################################################

    # overPlot <- overPlot +
    #   theme(axis.title.x = element_blank(),
    #         axis.text.x = element_blank(),
    #         axis.ticks.x = element_blank())
    # overPlot <- suppressMessages(
    #   plotly::ggplotly(overPlot)
    # )
  }

  overPlot

  # When moving to S4, maybe add the top 5-10 for a single file?
  # Add functionality to export a FASTA file of the sequences to the shiny app?
  # This will obviously work best under plotly as the names will be silly otherwise

}
