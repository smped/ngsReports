#' @title Plot GC Content for each file
#'
#' @description Plot the GC content for each file
#'
#' @details makes interactive heatmap for \code{GC content}
#'
#' For large datasets, subsetting by R1 or R2 reads may be helpful
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
#' @param counts \code{logical}. Display counts of GC content rather than frequency
#' @param pwfcols Object of class \code{\link{Pwfcol}} to give colours for pass, warning, and fail
#' values in plot
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param dendrogram \code{logical} redundant if \code{clusterNames} is \code{FALSE}
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param GCtheory \code{logical} default is \code{FALSE} to give the absolute value, set to \code{TRUE} to normalize
#' values of GC_Content by the theoretical values using \code{\link{gcTheoretical}}. \code{species} must be specified.
#' @param species \code{character} if \code{gcTheory} is \code{TRUE} its must be accompanied by a species
#' Currently supports Genome only (transcriptome to come). Species currently supported: A. lyrata, A. mellifera, A. thaliana,
#' B. taurus, C. elegans, C. familiaris, D. melanogaster, D. rerio, E. coli, G. aculeatus, G. gallus, H. sapiens, M. fascicularis,
#' M. furo, M. mulatta, M. musculus, O. sativa, P. troglodytes, R. norvegicus, S. cerevisiae, S scrofa, T. gondii,
#' T. guttata, V. vinifera. Use \code{ngsReports::genomes(ngsReports::gcTheoretical)} to display the corresponding names for
#' each species.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#'
#' @return A plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' barcodes <- c("ATTG", "CCGC", "CCGT", "GACC", "TTAT", "TTGG")
#' suffix <- c("R1_fastqc.zip", "R2_fastqc.zip")
#' fileList <- paste(rep(barcodes, each = 2), rep(suffix, times = 5), sep = "_")
#' fileList <- system.file("extdata", fileList, package = "fastqcReports")
#'
#' # Load the FASTQC data as a FastqcDataList
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot
#' plotGcHeatmap(fdl)
#'
#' # Using counts
#' plotGcHeatmap(fdl, counts = TRUE)
#'
#' @import tidyr
#' @import ggplot2
#' @import scales
#' @import plotly
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#'
#' @export plotGCHeatmapPlotly
plotGCHeatmapPlotly <- function(x, subset, counts = FALSE, pattern = "(.+)\\.(fastq|fq).*",
                          clusterNames = TRUE, pwfcols, GCtheory = FALSE, species = "Hsapiens"){
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if(GCtheory){
    spp <- ngsReports::genomes(ngsReports::gcTheoretical)
    if(!species %in% spp){
    stop(cat("Currently only supports the species", spp, sep = ", "))
      }
  }

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))

  pwfCols <- fastqcReports::pwf
  col <- getColours(pwfCols)

  x <- x[subset]
  df <- tryCatch(Per_sequence_GC_content(x))

  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>%
      dplyr::mutate(Freq = Count / sum(Count)) %>%
      dplyr::ungroup() %>%
      dplyr::select(Filename, GC_Content, Freq)

    # Add the observed mean at each value:
    mn <- dplyr::group_by(df, GC_Content) %>%
      dplyr::summarise(Freq = mean(Freq)) %>%
      dplyr::mutate(Freq = Freq / sum(Freq), # Normalise to 1 for a distribution
                    Filename = "Observed Mean")

    df <- dplyr::bind_rows(df, mn)
    df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

    if(GCtheory){
      df <- df %>% split(.["Filename"]) %>% lapply(function(x){
        gcTheoryDF <- ngsReports::getGC(ngsReports::gcTheoretical, species = Species, type = "Genome")
        x <- mutate(x, Freq = Freq - (gcTheoryDF$Genome/sum(gcTheoryDF$Genome)))
      }) %>% bind_rows()
    }


     if(clusterNames){
       df <- df %>% spread(GC_Content, Freq)
       rownames(df) <- as.matrix(df[1])
       df <- df  %>% as.data.frame()

       xx <- df %>% dplyr::select(-Filename)
       xx[is.na(xx)] <- 0
       clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
       row.ord <- order.dendrogram(clus)
       df <- df[row.ord,]
       df$Filename <- with(df, factor(Filename, levels=Filename, ordered=TRUE))
       dfLong <- df %>% tidyr::gather("GC_Content", "Frequency", 2:ncol(.))
       dfLong$GC_Content <- as.numeric(dfLong$GC_Content)
       dfLong$Frequency <- as.numeric(dfLong$Frequency)

       om_col <- data.frame(Status = NA, Category = NA, Filename = "Observed Mean")

       t <- getSummary(x) %>% dplyr::filter(Category == "Per sequence GC content")
       t <- bind_rows(om_col, t)
       t <- dplyr::full_join(df["Filename"], t, by = "Filename")
       t$Filename <- with(t, factor(Filename, levels=Filename))
       key <- t["Filename"]


       p <- ggplot2::ggplot(dfLong, aes(x = GC_Content, y = Filename)) + ggplot2::geom_tile(aes(fill = Frequency), color = "white", size = 30) +
         ggplot2::theme(panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank()) +
         ggplot2::scale_fill_gradientn(colours = viridis(50))

       d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
         ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                   panel.background = element_blank(),
                                                                   legend.position="none",
                                                                   axis.title=element_blank(),
                                                                   axis.text=element_blank(),
                                                                   axis.ticks=element_blank())
       d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


       GCheatmap <- subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "GC Content %"))
     }
if(!clusterNames){

  om_col <- data.frame(Status = NA, Category = NA, Filename = "Observed Mean")

  key <- df %>% spread(GC_Content, Freq) %>% select(Filename)
  t <- getSummary(x) %>% dplyr::filter(Category == "Per sequence GC content")
  t <- bind_rows(om_col, t)
  t <- dplyr::full_join(key, t, by = "Filename")
  t$Filename <- with(t, factor(Filename, levels=Filename))

  p <- ggplot2::ggplot(df, aes(x = GC_Content, y = Filename)) + ggplot2::geom_tile(aes(fill = Freq), color = "white", size = 30) +
    ggplot2::theme(panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank()) +
    ggplot2::scale_fill_gradientn(colours = viridis(50))

  d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
    ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(),
                                                              legend.position="none",
                                                              axis.title=element_blank(),
                                                              axis.text=element_blank(),
                                                              axis.ticks=element_blank())
  d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


  GCheatmap <- subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "GC Content %"))

}
  }
  else{

    # Add the observed mean at each value:
    mn <- dplyr::group_by(df, GC_Content) %>%
      dplyr::summarise(Count = mean(Count)) %>%
      dplyr::mutate(Filename = "Observed Mean")

    df <- dplyr::bind_rows(df, mn)
    df$Filename <- factor(df$Filename, levels = rev(unique(df$Filename)))

    if(GCtheory){
      df <- df %>% split(.["Filename"]) %>% lapply(function(x){
        gcTheoryDF <- ngsReports::getGC(ngsReports::gcTheoretical, species = Species, type = "Genome")
        x <- mutate(x, Count = Count - (gcTheoryDF$Genome/sum(gcTheoryDF$Genome)))
      }) %>% bind_rows()
    }

    if(clusterNames){
      df <- df %>% spread(GC_Content, Count)
      rownames(df) <- as.matrix(df[1])
      df <- df  %>% as.data.frame()

      xx <- df %>% dplyr::select(-Filename)
      xx[is.na(xx)] <- 0
      clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
      row.ord <- order.dendrogram(clus)
      df <- df[row.ord,]
      df$Filename <- with(df, factor(Filename, levels=Filename, ordered=TRUE))
      dfLong <- df %>% tidyr::gather("GC_Content", "Count", 2:ncol(.))
      dfLong$GC_Content <- as.numeric(dfLong$GC_Content)
      dfLong$Count <- as.numeric(dfLong$Count)

      om_col <- data.frame(Status = NA, Category = NA, Filename = "Observed Mean")

      t <- getSummary(x) %>% dplyr::filter(Category == "Per sequence GC content")
      t <- bind_rows(om_col, t)
      t <- dplyr::full_join(df["Filename"], t, by = "Filename")
      t$Filename <- with(t, factor(Filename, levels=Filename))
      key <- t["Filename"]


      p <- ggplot2::ggplot(dfLong, aes(x = GC_Content, y = Filename)) + ggplot2::geom_tile(aes(fill = Count), color = "white", size = 30) +
        ggplot2::theme(panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank()) +
        ggplot2::scale_fill_gradientn(colours = viridis(50))

      d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
        ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(),
                                                                  legend.position="none",
                                                                  axis.title=element_blank(),
                                                                  axis.text=element_blank(),
                                                                  axis.ticks=element_blank())
      d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


      GCheatmap <- subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "GC Content %"))
    }
    if(!clusterNames){

      om_col <- data.frame(Status = NA, Category = NA, Filename = "Observed Mean")

      key <- df %>% spread(GC_Content, Count) %>% select(Filename)
      t <- getSummary(x) %>% dplyr::filter(Category == "Per sequence GC content")
      t <- bind_rows(om_col, t)
      t <- dplyr::full_join(key, t, by = "Filename")
      t$Filename <- with(t, factor(Filename, levels=Filename))

      p <- ggplot2::ggplot(df, aes(x = GC_Content, y = Filename)) + ggplot2::geom_tile(aes(fill = Count), color = "white", size = 30) +
        ggplot2::theme(panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank()) +
        ggplot2::scale_fill_gradientn(colours = viridis(50))

      d <- ggplot2::ggplot(t, aes(x = 1, y = Filename, key = key)) + ggplot2::geom_tile(aes(fill = Status)) +
        ggplot2::scale_fill_manual(values = col) + ggplot2::theme(panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(),
                                                                  legend.position="none",
                                                                  axis.title=element_blank(),
                                                                  axis.text=element_blank(),
                                                                  axis.ticks=element_blank())
      d <- plotly::ggplotly(d, tooltip = c("Status", "Filename"))


      GCheatmap <- subplot(d, p, widths = c( 0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "GC Content %"))

    }

  }
 GCheatmap
}
