#' @title Plot the per base content as a heatmap
#'
#' @description Plot the Per Base content for a set of FASTQC files.
#' Informative plot where per base sequence content (%A, %T, %G, %C),
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly \code{logical}. Generate an interactive plot using plotly
#' @param plotType \code{character}. Type of plot to generate. Must be "line" or "heatmap"
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
#' plotSequenceContent(fdl)
#'
#'
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom scales rescale
#' @importFrom grDevices rgb
#' @importFrom magrittr %>%
#'
#' @export

plotSequenceContent <- function(x, usePlotly = FALSE, labels, plotType = "heatmap"){

  stopifnot(plotType %in% c("heatmap", "line"))
  # Get the SequenceContent
  df <- tryCatch(Per_base_sequence_content(x))
  df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))

  # Drop the suffix, or check the alternate labels
  if (missing(labels)){
    labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
                        names = fileName(x))
  }
  else{
    if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
  }
  if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")

  if(plotType == "heatmap"){

  maxBase <- max(vapply(c("A", "C", "G", "T"), function(x){max(df[[x]])}, numeric(1)))
  df$colour <- with(df, rgb(floor(`T`) / maxBase, floor(A) / maxBase, floor(C) / maxBase, 1 - floor(G) / maxBase))

  basicStat <- Basic_Statistics(x)[c("Filename", "Longest_sequence")]

  df <- dplyr::right_join(df, basicStat, by = "Filename")
  df <- df[c("Filename", "Start", "colour", "Longest_sequence", "A", "C", "G", "T")]

  df <- split(df, f = df$Filename) %>%
    lapply(function(x){
      dfFill <- data.frame(Start = 1:x[["Longest_sequence"]][1])
      x <- dplyr::right_join(x, dfFill, by = "Start") %>%
        zoo::na.locf()
    }) %>%
    dplyr::bind_rows()
  df$Start <- as.integer(df$Start)
  tileCols <- unique(df$colour)
  names(tileCols) <- unique(df$colour)

  sequenceContentHeatmap <- ggplot(df,
                                   aes_string(x = "Start", y = "Filename",
                                       fill = "colour", A = "A", C = "C", G = "G", T = "T")) +
    geom_tile() +
    scale_fill_manual(values = tileCols) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), labels = labels) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "black"))

  if (usePlotly){
    suppressMessages(plotly::ggplotly(sequenceContentHeatmap, tooltip = c("A", "C", "G", "T", "Start", "Filename")))
  }
  else{
    sequenceContentHeatmap
  }}
  else{
    df$Filename <- labels[df$Filename]
    df2 <- df[!colnames(df) == "Base"]
    df2 <- melt(df2, id.vars = c("Filename", "Start"))
    colnames(df2) <- c("Filename", "Start", "Base", "Percent")
    df2$Base <- factor(df2$Base, levels = unique(df2$Base))

    #set colours
    baseCols <- factor(c(`T`="red", G = "black", A = "green", C = "blue"))

    sequenceContentHeatmap <- ggplot(df2, aes_string(x = "Start", y = "Percent", colour = "Base")) +
      geom_line() +
      facet_wrap(~Filename) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      guides(fill = FALSE) +
      labs(x = "Position in read (bp)",
           y = "Percent (%)") +
      theme_bw() +
      scale_colour_manual(values = baseCols)
    if(usePlotly){
      ggplotly(sequenceContentHeatmap)
    }else{
      sequenceContentHeatmap
    }
  }
}

