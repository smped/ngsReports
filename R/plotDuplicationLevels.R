#' @title Plot the combined Sequence_Duplication_Levels information
#'
#' @description Plot the Sequence_Duplication_Levels information for a set of FASTQC reports
#'
#' @details
#' This extracts the Sequence_Duplication_Levels from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileNames
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param type A regular expression used to filter which value is plotted.
#' Patterns should match one of \code{\% Total sequences} or \code{\% Deduplicated sequences}
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
#' # Draw the default plot
#' plotDuplicationLevels(fdl)
#'
#' # Plot the % Total Sequences for the R1 files only
#' r1 <- grepl("R1", fileNames(fdl))
#' plotDuplicationLevels(fdl, subset = r1, type = "Total")
#'
#' 
#' @importFrom stringr str_detect
#' @importFrom stringr str_to_title
#' @importFrom reshape2 melt
#' @importFrom dplyr filter
#'
#' @export
plotDuplicationLevels <- function(x, subset, type = ".+",
                                  trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*"){

  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(is.character(type))

  x <- x[subset]
  df <- tryCatch(Sequence_Duplication_Levels(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

    # Melt for a convenient long form
  df <- reshape2::melt(df, id.vars = c("Filename", "Duplication_Level"),
                       variable.name = "Type", value.name = "Percent")

  # Ensure the Duplication_level plots in the correct order
  df$Duplication_Level <- factor(df$Duplication_Level,
                                 levels = c(1:9, paste0(">", c(10, 50, 100, 500, "1k", "5k", "10k+"))))

  # Tidy the output to look like the original report
  df$Type <- gsub("Percentage_of_", "% ", df$Type)
  df$Type <- stringr::str_to_title(df$Type)
  df$Type <- paste(df$Type, "sequences")

  # Restrict to a given type if requested
  df <- dplyr::filter(df, grepl(type, Type))
  if(nrow(df) == 0) stop("Invalid type selected")

  dupPlot <- ggplot(df, aes(x = as.integer(Duplication_Level), y = Percent)) +
    annotate("rect",
                      xmin = seq(1.5, 17, by = 2), xmax = seq(2.5, 17, by = 2),
                      ymin = 0, ymax = Inf, fill = "grey80", alpha = 0.5) +
    geom_line(aes(colour = Filename, group = Filename)) +
    facet_wrap(~Type, ncol = 1) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                                breaks = seq(0, 100, by = 20)) +
    scale_x_continuous(breaks = seq_along(levels(df$Duplication_Level)),
                                labels = levels(df$Duplication_Level),
                                limits = c(0, length(levels(df$Duplication_Level))) + 0.5,
                                expand = c(0, 0)) +
    labs(x = "Duplication Level", y = "Percent (%)") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())

    # And draw the plot
    dupPlot
}
