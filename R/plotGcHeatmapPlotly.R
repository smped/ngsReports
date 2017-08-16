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
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileName.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param clusterNames \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using heirachial clustering
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param GCtheory \code{logical} default is \code{FALSE} to give the true GC content%, set to \code{TRUE} to normalize
#' values of GC_Content by the theoretical values using \code{\link{gcTheoretical}}. \code{species} must be specified.
#' @param GCtheoryType \code{"character"} Select type of data to normalize GC content% agianst accepts either "Genome" or
#' "Transcriptome" Default is "Genome"
#' @param species \code{character} if \code{gcTheory} is \code{TRUE} its must be accompanied by a species
#' Currently supports Genome only (transcriptome to come). Species currently supported: A. lyrata, A. mellifera, A. thaliana,
#' B. taurus, C. elegans, C. familiaris, D. melanogaster, D. rerio, E. coli, G. aculeatus, G. gallus, H. sapiens, M. fascicularis,
#' M. furo, M. mulatta, M. musculus, O. sativa, P. troglodytes, R. norvegicus, S. cerevisiae, S scrofa, T. gondii,
#' T. guttata, V. vinifera. Use \code{ngsReports::genomes(ngsReports::gcTheoretical)} to display the corresponding names for
#' each species.
#' @param trimNames \code{logical}. Capture the text specified in \code{pattern} from fileName
#' @param dendrogram \code{logical} redundant if \code{clusterNames} and \code{usePlotly} are \code{FALSE}.
#' if both \code{clusterNames} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#'
#' @return A plotly object
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
#' plotGcHeatmap(fdl)
#'
#' # Using counts
#' plotGcHeatmap(fdl, counts = TRUE)
#'
#'
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#'
#' @export plotGCHeatmapPlotly
plotGCHeatmapPlotly <- function(x, subset, counts = FALSE, pattern = "(.+)\\.(fastq|fq).*",
                                clusterNames = FALSE, pwfCols,
                                GCtheory = FALSE, GCtheoryType = "Genome", species = "Hsapiens",
                                trimNames = TRUE, usePlotly = FALSE,
                                dendrogram = FALSE){
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if(GCtheory & GCtheoryType == "Genome"){
    spp <- ngsReports::genomes(ngsReports::gcTheoretical)
    if(!species %in% spp$Name){
      stop(cat("Currently only supports genomes for", spp$Name, sep = ", "))
    }
  }

  if(GCtheory & GCtheoryType == "Transcriptome"){
    spp <- ngsReports::transcriptomes(ngsReports::gcTheoretical)
    if(!species %in% spp$Name){
      stop(cat("Currently only supports transcriptomes for", spp$Name, sep = ", "))
    }
  }

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))

  if (missing(pwfCols)) pwfCols <- ngsReports::pwf
  col <- ngsReports::getColours(pwfCols)

  x <- x[subset]
  df <- tryCatch(Per_sequence_GC_content(x))

  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }else{
    pattern <- ""
  }


  if (!counts){

    # Summarise to frequencies & initialise the plot
    df <- dplyr::group_by(df, Filename) %>%
      dplyr::mutate(Freq = Count / sum(Count)) %>%
      dplyr::ungroup() %>%
      dplyr::select(Filename, GC_Content, Value = Freq)

    if(GCtheory){
      df <- df %>% split(.["Filename"]) %>% lapply(function(x){
        gcTheoryDF <- ngsReports::getGC(ngsReports::gcTheoretical, name = species, type = GCtheoryType)
        x <- dplyr::mutate(x, Value = abs(Value - unlist(gcTheoryDF[species])))
      }) %>% dplyr::bind_rows(.)
    }
  }else{
    df <- dplyr::select(df, Filename, GC_Content, Value = Count)
  }

  if(clusterNames){
    df <- reshape2::dcast(df, Filename ~ GC_Content)
    xx <- dplyr::select(df, -Filename)
    xx[is.na(xx)] <- 0
    clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
    row.ord <- order.dendrogram(clus)
    df <- df[row.ord,]
    df <- reshape2::melt(df, id.vars = "Filename", variable.name = "GC_Content", value.name = "Value")
  }

  if(!counts){
    df <- dplyr::mutate(df, Frequency = as.numeric(Value),
                        GC_Content = as.integer(GC_Content),
                        Filename = factor(Filename, levels = unique(Filename)))
    GCheatmap <- ggplot(df, aes(x = GC_Content, y = Filename, fill = Frequency))
  }else{
    df <- dplyr::mutate(df, Counts = as.numeric(Value),
                        GC_Content = as.integer(GC_Content),
                        Filename = factor(Filename, levels = unique(Filename)))
    GCheatmap <- ggplot(df, aes(x = GC_Content, y = Filename, fill = Counts))
  }

  GCheatmap <- GCheatmap + geom_tile() +
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_fill_gradientn(colours = viridisLite::inferno(50))

  if(usePlotly){

    GCheatmap <- GCheatmap + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())


    t <- getSummary(x) %>% dplyr::filter(Category == "Per sequence GC content")
    t <- dplyr::mutate(t, FilenameFull = Filename,
                       Filename = gsub(pattern[1], "\\1", t$Filename),
                       Filename = factor(Filename, levels = unique(df$Filename)))
    t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
    key <- t$FilenameFull

    sideBar <- ggplot(t, aes(x = 1, y = Filename, key = key)) + geom_tile(aes(fill = Status)) +
      scale_fill_manual(values = col) + theme(panel.grid.minor = element_blank(),
                                              panel.background = element_blank(),
                                              legend.position="none",
                                              axis.title=element_blank(),
                                              axis.text=element_blank(),
                                              axis.ticks=element_blank())
    sideBar <- plotly::ggplotly(sideBar, tooltip = c("Status", "Filename"))

    #plot dendrogram
    if(dendrogram){
      ggdend <- function(df) {
        ggplot() +
          geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) + ggdendro::theme_dendro()
      }

      dx <- ggdendro::dendro_data(clus)
      dendro <- ggdend(dx$segments) + coord_flip() +
        scale_y_reverse(expand = c(0, 1)) + scale_x_continuous(expand = c(0,1))


      GCheatmap <- plotly::subplot(dendro, sideBar, GCheatmap, widths = c(0.3, 0.1,0.6), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis3 = list(title = "GC Content (%)"))
    }else{
      GCheatmap <- plotly::subplot(sideBar, GCheatmap, widths = c(0.1,0.9), margin = 0, shareY = TRUE) %>% plotly::layout(xaxis2 = list(title = "GC Content (%)"))
    }

  }else{
    GCheatmap
  }
  GCheatmap
}
