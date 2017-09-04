#' @title Plot the Per Sequence GC Content
#'
#' @description Plot the Per Sequence GC Content for a set of FASTQC files
#'
#' @details
#' Draws a plot of all GC Content overlaid across a single plot.
#' The mean from the observed data at each value is overlaid as a dotted black line.
#' Addition of a theoretical distribution has not yet been implemented.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param counts \code{logical}. Plot the counts from each file if \code{counts = TRUE}.
#' If \code{counts = FALSE} the frequencies will be plotted
#' @param labels An optional named factor of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#' @return A plotly object
#'
#' @return A ggplot2 object
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
#' # The default plot
#' plotGcContent(fdl)
#'
#' # Plot the R1 files using counts
#' r1 <- grepl("R1", fileName(fdl))
#'
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom plotly ggplotly
#'
#' @export
plotGcContent <- function(x, usePlotly = FALSE, labels,
                          GCtheory = FALSE, GCtheoryType = "Genome",
                          species = "Hsapiens", GCobject, counts = FALSE, ...){

  df <- tryCatch(Per_sequence_GC_content(x))

  # Get the correct y-axis label
  ylab <- c("Frequency", "Count")[counts + 1]

  # Remove zero counts
  # df <- dplyr::filter(df, Count > 0)

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  # Get any arguments for dotArgs that have been set manually
  dotArgs <- list(...)
  allowed <- names(formals(ggplot2::theme))
  keepArgs <- which(names(dotArgs) %in% allowed)
  userTheme <- c()
  if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

  # Tidy up the GC content variables
  if(GCtheory & missing(GCobject)){
    GCobject <- ngsReports::gcTheoretical
  }

  if(GCtheory & GCtheoryType == "Genome"){
    spp <- genomes(GCobject)
    if(!species %in% spp$Name){
      stop(cat("Currently only supports genomes for", spp$Name, sep = ", "))
    }
  }

  if(GCtheory & GCtheoryType == "Transcriptome"){
    spp <- transcriptomes(GCobject)
    if(!species %in% spp$Name){
      stop(cat("Currently only supports transcriptomes for", spp$Name, sep = ", "))
    }
  }


  df$Filename <- labels[df$Filename]


  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>%
      dplyr::mutate(Freq = Count / sum(Count)) %>%
      dplyr::ungroup() %>%
      dplyr::select(Filename, GC_Content, Freq)

     if(GCtheory){
        gcTheoryDF <- getGC(GCobject, name = species, type = GCtheoryType)
        names(gcTheoryDF)[names(gcTheoryDF) == species] <- "Freq"

        gcTheoryDF$Filename <- paste("Theoretical GC for", species)

        df <- bind_rows(df, gcTheoryDF)

     }
    df$Filename <- factor(df$Filename, levels = unique(df$Filename))
    gcPlot <- ggplot(df, aes(x = GC_Content, y = Freq, colour = Filename)) +
      geom_line()

  }
  else{

    # Initialise the plot using counts
    gcPlot <- ggplot(df, aes(x = GC_Content, y = Count, colour = Filename)) +
      geom_line()

  }

  # Add the rest of the plotting detail
  gcPlot <- gcPlot +
    scale_colour_discrete(labels = labels) +
    ylab(ylab) +
    theme_bw()
  if (!is.null(userTheme)) gcPlot <- gcPlot + userTheme

  if(usePlotly){
    if(!counts) value <- "Freq"
    else value <- "Count"
    gcPlot <- gcPlot + labs(colour = "")
    gcPlot <- ggplotly(gcPlot, tooltip = c("GC_Content", value))
  }

  # Draw the plot
  gcPlot
}
