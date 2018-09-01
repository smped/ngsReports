#' @title Plot Overrepresented Kmers
#'
#' @description Plot Overrepresented Kmers
#'
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param n \code{numeric}. The number of Kmers to show.
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default.
#' @param usePlotly \code{logical} Default \code{FALSE} will render using ggplot.
#' If \code{TRUE} plot will be rendered with plotly
#' @param ... Used to pass various potting parameters to theme.
#' Can also be used to set size and colour for box outlines.
#' @param lineWidth Passed to \code{geom_line(size = lineWidth)}
#' @param pal The colour palette.
#' If the vector supplied is less than n, \code{grDevices::colorRampPalette()} will be used
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param heatCol Colour palette used for the heatmap. Default is \code{inferno} from the package \code{viridris}
#'
#' @return A standard ggplot2 object or an interactive plotly object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#' plotKmers(fdl[[1]])
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr desc
#' @import ggplot2
#'
#' @name plotKmers
#' @rdname plotKmers-methods
#' @export
setGeneric("plotKmers",function(x, usePlotly = FALSE, ...){standardGeneric("plotKmers")})
#' @aliases plotKmers,character
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "character",
          function(x, usePlotly = FALSE, ...){
              x <- getFastqcData(x)
              plotKmers(x, usePlotly,...)
          }
)
#' @aliases plotKmers,FastqcFile
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "FastqcFile",
          function(x, usePlotly = FALSE, ...){
              x <- getFastqcData(x)
              plotKmers(x, usePlotly,...)
          }
)
#' @aliases plotKmers,FastqcFileList
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, ...){
              x <- getFastqcData(x)
              plotKmers(x, usePlotly,...)
          }
)
#' @aliases plotKmers,FastqcData
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "FastqcData",
          function(x, usePlotly = FALSE,  n = 6, labels, ..., lineWidth = 0.5,
                   pal = c("red", "blue", "green", "black", "magenta", "yellow")){
              
              # Get the basic data frame
              df <- Kmer_Content(x)
              
              if (!length(df)) {
                  #stop("No Overrep kmers")
                  kMerPlot <- emptyPlot("No Adapter Content Module Detected")
                  
                  if(usePlotly) kMerPlot <- ggplotly(kMerPlot, tooltip = "")
                  return(kMerPlot)
              }
              
              
              # Drop the suffix, or check the alternate labels
              if (missing(labels)){
                  labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
              }
              else{
                  if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
              }
              if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
              df$Filename <- labels[df$Filename]
              
              # Get the top kMers
              o <- order(df$PValue, 1/df$Count)
              allK <- unique(df$Sequence)
              n <- tryCatch(as.integer(n))
              n <- min(length(allK), n)
              topK <- unique(df$Sequence[o])[seq_len(n)]
              # Tidy up the data
              df <- dplyr::filter(df, Sequence %in% topK)
              colnames(df) <- gsub("Max_Obs/Exp_Position", "Base", colnames(df))
              colnames(df) <- gsub("Obs/Exp_Max", "Value", colnames(df))
              df$Position <- gsub("([0-9]*)-[0-9]*", "\\1", df$Base)
              df$Position <- as.integer(df$Position)
              df <- df[c("Filename", "Sequence", "Value", "Position")]
              
              # Set the x-axis for plotting
              # As there will be significant gaps in this data,
              # the bins used need to be obtained from another slot.
              # The most complete will be Per_base_sequence_quality
              # These values can then be incorporated in the final df for accurate plotting & labelling
              refForX <- unique(Per_base_sequence_quality(x)$Base)
              refForX <- tibble::tibble(Base = as.character(refForX),
                                        Position = gsub("([0-9]*)-[0-9]*", "\\1", Base))
              refForX$Position <- as.integer(refForX$Position)
              
              # In order to get a line plot, zero values need to be added to the missing positions
              # The above reference scale for X will be used to label the missing values
              # Include all files to ensure all appear in the final plot
              zeros <- expand.grid(list(Filename = df$Filename,
                                        Sequence = topK,
                                        Value = 0,
                                        Position = refForX$Position),
                                   stringsAsFactors = FALSE)
              
              # After the bind_rows, duplicate values will exist at some positions
              # Spuriously introduced zeros need to be removed
              df <- dplyr::bind_rows(df, zeros)
              df <- dplyr::arrange(df, Filename, Sequence, Position, desc(Value))
              df <- dplyr::distinct(df, Filename, Sequence, Position, .keep_all = TRUE)
              df <- dplyr::left_join(df, refForX, by = "Position")
              
              # Set the Sequence as a factor based on the first position it appears
              # This way the colours will appear in order in the guide as well as the plot
              kMerLevels <- dplyr::arrange(df, Position, Sequence) %>%
                  dplyr::filter(Value > 0) %>%
                  magrittr::extract2("Sequence") %>%
                  unique()
              df$Sequence <- factor(df$Sequence, levels = kMerLevels)
              df <- droplevels(df)
              
              # Set the plotting params
              yMax <- max(df$Value)*1.05
              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              # And the colours
              if (n < length(pal)) pal <- pal[seq_len(n)]
              if (n > length(pal)) pal <- grDevices::colorRampPalette(pal)(n)
              
              # Now draw the basic plots
              kMerPlot <- ggplot(df,
                                 aes_string(x =" Position", y = "Value", colour = "Sequence")) +
                  geom_line(size = lineWidth) +
                  facet_wrap(~Filename) +
                  scale_x_continuous(breaks = refForX$Position,
                                     labels = refForX$Base,
                                     expand = c(0.02, 0)) +
                  scale_y_continuous(limits = c(0, yMax), expand = c(0, 0)) +
                  scale_colour_manual(values = pal) +
                  theme_bw() +
                  labs(x = "Position in read (bp)",
                       y = expression(paste(log[2], " Obs/Exp")),
                       colour = c())
              
              # Check for binned x-axis values to decied whether to rotate x-axis labels
              # This should be clear if there are more than 2 characters in the plotted labels
              binned <- any(grepl("-", df$Base))
              if (binned) kMerPlot <- kMerPlot + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))
              if (!is.null(userTheme)) kMerPlot <- kMerPlot + userTheme
              
              if(usePlotly){
                  kMerPlot <- kMerPlot +
                      ylab("Log2 Obs/Exp") +
                      xlab("")
                  kMerPlot <- suppressMessages(
                      plotly::subplot(plotly::plotly_empty(), kMerPlot, widths = c(0.14,0.86)) %>% 
                          layout(yaxis2 = list(title = "Log(2) Obs/Exp"),
                                 legend = list(orientation = "v", x = 1, y = 1)))
              }
              kMerPlot
          }
)
#' @aliases plotKmers,FastqcDataList
#' @rdname plotKmers-methods
#' @export
setMethod("plotKmers", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, cluster = FALSE,
                   dendrogram = FALSE, pwfCols, heatCol = inferno(50), ...){
              
              df <- Kmer_Content(x)
              
              if (!length(df)) {
                  #stop("No Overrep kmers")
                  kMerPlot <- emptyPlot("No Overrepresented Kmers Detected")
                  
                  if(usePlotly) kMerPlot <- ggplotly(kMerPlot, tooltip = "")
                  return(kMerPlot)
              }
              
              # Sort out the colours
              if (missing(pwfCols)) pwfCols <- ngsReports::pwf
              stopifnot(isValidPwf(pwfCols))
              
              # Drop the suffix, or check the alternate labels
              if (missing(labels)){
                  labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)), names = fileName(x))
              }
              else{
                  if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
              }
              if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated Totals")
              
              
              colnames(df) <- gsub("Max_Obs/Exp_Position", "Base", colnames(df))
              colnames(df) <- gsub("Obs/Exp_Max", "Total", colnames(df))
              df$Position <- gsub("([0-9]*)-[0-9]*", "\\1", df$Base)
              df$Position <- as.integer(df$Position)
              
              
              # Now simply add the kmers at each position. Not superinformative...
              df <- dplyr::group_by(df, Filename, Position)
              df <- dplyr::summarise(df, Total = sum(Count))
              df <- dplyr::ungroup(df)
              
              # Setup the x-axis
              refForX <- unique(Per_base_sequence_quality(x)$Base)
              refForX <- tibble::tibble(Base = as.character(refForX),
                                        Position = gsub("([0-9]*)-[0-9]*", "\\1", Base))
              refForX$Position <- as.integer(refForX$Position)
              
              # Setup a grid of NAs
              nas <- expand.grid(list(Filename = unique(df$Filename),
                                      Total = NA,
                                      Position = unique(refForX$Position)),
                                 stringsAsFactors = FALSE)
              
              df <- dplyr::bind_rows(df, nas)
              df <- dplyr::arrange(df, Filename, Position, desc(Total))
              df <- dplyr::distinct(df, Filename, Position, .keep_all = TRUE)
              
              
              # casting before going into the loop makes getting th file names in order much easier
              df <- reshape2::dcast(df, Filename ~ Position, value.var = "Total")
              
              #cluster
              if(cluster){
                  xx <- df[!colnames(df) == "Filename"]
                  xx[is.na(xx)] <- 0
                  clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
                  row.ord <- order.dendrogram(clus)
                  df <- df[row.ord,]
              }
              
              key <- df$Filename
              df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Position", value.name = "Total")
              # I dont like this
              df$Position <- as.integer(as.character(df$Position))
              # Set Filenames as factor, make key for sidebar and change Filenames to match provided/
              # default labels
              df$Filename <- labels[df$Filename]
              df$Filename <- factor(df$Filename, levels = unique(df$Filename))
              
              # join x axis ticks and df
              df <- dplyr::left_join(df, refForX, by = "Position")
              
              # Set up for geom_tile
              df$End <- gsub("[0-9]*-([0-9]*)", "\\1", df$Base)
              df$End <- as.integer(df$End)
              df$`Middle of Bin` <- (df$Position + df$End) / 2
              df$Width <- df$End - df$Position
              df$Width[df$Width == 0] <- 1
              
              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              kMerPlot <- ggplot(df, aes_string(x = "`Middle of Bin`", y = "Filename", fill = "Total", width = "Width")) +
                  geom_tile() +
                  scale_x_continuous(expand = c(0.02, 0)) +
                  scale_fill_gradientn(colors = heatCol, na.value = "white") +
                  theme_bw()
              
              if (!is.null(userTheme)) kMerPlot <- kMerPlot + userTheme
              
              if (usePlotly){
                  kMerPlot <- kMerPlot + theme(axis.text.y = element_blank())
                  t <- getSummary(x)
                  t <- t[t$Category == "Kmer Content",]
                  t$Filename <- labels[t$Filename]
                  t$Filename <- factor(t$Filename, levels = levels(df$Filename))
                  
                  # Key was not in same order as t
                  t <- t[order(t$Filename, levels(t$Filename)),]
                  
                  
                  sideBar <- makeSidebar(status = t, key = key, pwfCols = pwfCols)
                  
                  if (cluster && dendrogram){
                      dx <- ggdendro::dendro_data(clus)
                      dendro <- ggdend(dx$segments) +
                          coord_flip() +
                          scale_y_reverse(expand = c(0, 0)) +
                          scale_x_continuous(expand = c(0, 0.5))
                      
                      kMerPlot <- suppressWarnings(
                          suppressMessages(
                              plotly::subplot(dendro, sideBar, kMerPlot, widths = c(0.08, 0.09,0.83),
                                              margin = 0.001, shareY = TRUE, shareX = TRUE)
                          ))
                  }
                  else{
                      # Return the plot
                      kMerPlot <- suppressWarnings(
                          suppressMessages(
                              plotly::subplot(plotly::plotly_empty(),
                                              sideBar, kMerPlot, widths = c(0.08, 0.09, 0.83), margin = 0.001, shareY = TRUE) %>%
                                  plotly::layout(annotations = list(text = "Filename", showarrow = FALSE,
                                                                    textangle = -90))
                          ))
                  }
                  
                  kMerPlot <- kMerPlot %>%
                      plotly::layout(xaxis3 = list(title = "Position in Read (bp)"),
                                     margin = list(b = 50))
              }
              
              kMerPlot
              
          }
)
