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
#' @param pwfCols Object of class \code{\link{PwfCols}} to give colours for pass, warning, and fail
#' values in plot
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param dendrogram \code{logical} redundant if \code{cluster} is \code{FALSE}
#' if both \code{cluster} and \code{dendrogram} are specified as \code{TRUE} then the dendrogram
#' will be displayed.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A ggplot2 object
#'
#' @examples
#'
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # The default plot
#' plotSequenceContent(fdl)
#'
#'
#' @importFrom scales rescale
#' @importFrom grDevices rgb
#' @importFrom dplyr mutate_at vars funs
#' @importFrom tidyselect one_of
#' @import ggplot2
#'
#' @name plotSequenceContent
#' @rdname plotSequenceContent-methods
#' @export
setGeneric("plotSequenceContent",function(x, usePlotly = FALSE, labels, ...){standardGeneric("plotSequenceContent")})
#' @aliases plotSequenceContent,character
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "character",
          function(x, usePlotly = FALSE, labels, ...){
              x <- getFastqcData(x)
              plotSequenceContent(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceContent,FastqcFile
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcFile",
          function(x, usePlotly = FALSE, labels, ...){
              x <- getFastqcData(x)
              plotSequenceContent(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceContent,FastqcFileList
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcFileList",
          function(x, usePlotly = FALSE, labels, ...){
              x <- getFastqcData(x)
              plotSequenceContent(x, usePlotly, labels, ...)
          }
)
#' @aliases plotSequenceContent,FastqcData
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcData",
          function(x, usePlotly = FALSE, labels, ...){
              
              # Get the SequenceContent
              df <- Per_base_sequence_content(x)
              
              if (!length(df)) {
                  scPlot <- emptyPlot("No Sequence Content Module Detected")
                  if(usePlotly) scPlot <- ggplotly(scPlot, tooltip = "")
                  return(scPlot)
              }
              
              df$Position <- gsub("([0-9]*)-[0-9]*", "\\1", df$Base)
              df$Position <- as.integer(df$Position)
              
              # Drop the suffix, or check the alternate labels
              labels <- setLabels(df, labels, ...)
              
              df$Filename <- labels[df$Filename]
              df <- df[!colnames(df) == "Base"]
              df <- tidyr::gather(df, key = "Base", value = "Percent", 
                            tidyselect::one_of(c("G", "A", "T", "C")))
              df$Base <- factor(df$Base, levels = c("T", "C", "A", "G"))
              df$Percent <- round(df$Percent, 2)
              
              #set colours
              baseCols <- c(`T`= "red", G = "black", A = "green", C = "blue")
              
              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              add_percent <- function(x){paste0(x, "%")}
              xLab <- "Position in read (bp)"
              yLab <- "Percent"
              scPlot <- ggplot(df, aes_string(x = "Position", y = "Percent", colour = "Base")) +
                  geom_line() +
                  facet_wrap(~Filename) +
                  scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                                     labels = add_percent) +
                  scale_x_continuous(expand = c(0, 0)) +
                  guides(fill = FALSE) +
                  labs(x = xLab, y = yLab) +
                  theme_bw() +
                  scale_colour_manual(values = baseCols)
              
              if(usePlotly){
                  
                  scPlot <- suppressMessages(
                      suppressWarnings(
                          plotly::subplot(plotly::plotly_empty(), scPlot, 
                                          widths = c(0.14,0.86))
                      )
                  )
                  
                  scPlot <- plotly::layout(scPlot,
                                           xaxis2 = list(title = xLab), 
                                           yaxis2 = list(title = yLab))
              }
              
              scPlot
              
          }
)
#' @aliases plotSequenceContent,FastqcDataList
#' @rdname plotSequenceContent-methods
#' @export
setMethod("plotSequenceContent", signature = "FastqcDataList",
          function(x, usePlotly = FALSE, labels, plotType = c("heatmap", "line"), pwfCols,
                   cluster = FALSE, dendrogram = TRUE, ...){
              
              # Get the SequenceContent
              df <- Per_base_sequence_content(x)
              
              if (!length(df)) {
                  scPlot <- emptyPlot("No Sequence Content Module Detected")
                  if(usePlotly) scPlot <- ggplotly(scPlot, tooltip = "")
                  return(scPlot)
              }
              
              df$Start <- gsub("([0-9]*)-[0-9]*", "\\1", df$Base)
              df$End <- gsub("[0-9]*-([0-9]*)", "\\1", df$Base)
              df <- dplyr::mutate_at(df, vars(Start, End), funs(as.integer))
              
              plotType <- match.arg(plotType)
              if (missing(pwfCols)) pwfCols <- ngsReports::pwf
              
              # Drop the suffix, or check the alternate labels
              labels <- setLabels(df, labels, ...)
              
              # Get any arguments for dotArgs that have been set manually
              dotArgs <- list(...)
              allowed <- names(formals(ggplot2::theme))
              keepArgs <- which(names(dotArgs) %in% allowed)
              userTheme <- c()
              if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
              
              if (plotType == "heatmap"){
                  
                  # Round to 2 digits to reduce the complexity of the colour palette
                  df <- dplyr::mutate_at(df, vars(one_of(c("A", "C", "G", "T"))),
                                                  funs(round), digits = 2)
                  maxBase <- max(vapply(c("A", "C", "G", "T"), 
                                        function(x){max(df[[x]])}, numeric(1)))
                  # Set the colours, using opacity for G
                  df$opacity <- 1 - df$G / maxBase
                  df$colour <- with(df, rgb(red = `T` * opacity / maxBase,
                                            green = A * opacity / maxBase,
                                            blue = C * opacity / maxBase))
                  
                  basicStat <- Basic_Statistics(x)[c("Filename", "Longest_sequence")]
                  df <- dplyr::right_join(df, basicStat, by = "Filename")
                  df <- df[c("Filename", "Start", "End", "colour", 
                             "Longest_sequence", "A", "C", "G", "T")]
                  # df$Position <- as.integer(df$Start)
                  # df$Filename <- labels[df$Filename]
                  
                  if (dendrogram && !cluster){
                      message("cluster will be set to TRUE when dendrogram = TRUE")
                      cluster <- TRUE
                  }
                  
                  # Now define the order for a dendrogram if required
                  key <- names(labels)
                  if (cluster){
                      df_gath <- tidyr::gather(df, key = "Base", value = "Percent", 
                                               one_of(c("A", "C", "G", "T")))
                      df_gath$Start <- paste(df_gath$Start, df_gath$Base, sep = "_")
                      cols <- c("Filename", "Start", "Percent")
                      clusterDend <- setClusters(df = df_gath[cols], 
                                                 rowVal = "Filename", 
                                                 colVal = "Start", 
                                                 value = "Percent")
                      key <- labels(clusterDend)
                  }
                  # Now set everything as factors
                  df$Filename <- factor(labels[df$Filename], 
                                        levels = labels[key])
                  # Define the colours as named colours (name = colour)
                  tileCols <- unique(df$colour)
                  names(tileCols) <- unique(df$colour)
                  # Define the tile locations
                  df$y <- as.integer(df$Filename)
                  df$ymax <- as.integer(df$Filename) + 0.5
                  df$ymin <- df$ymax - 1
                  df$xmax <- df$End + 0.5
                  df$xmin <- df$Start - 1
                  df$Window <- paste(df$Start, "-", df$End, "bp")
                  
                  xLab <- "Position in read (bp)"
                  yLab <- "Filename"
                  
                  scPlot <- ggplot(df, 
                                   aes_string(y = "y",
                                              fill = "colour", 
                                              A = "A", C = "C", G = "G", `T` = "T", 
                                              Filename = "Filename", Window = "Window")) + 
                      geom_rect(aes_string(xmin = "xmin", xmax = "xmax", 
                                           ymin = "ymin", ymax = "ymax"),
                                colour = "#00000000") +
                      scale_fill_manual(values = tileCols) +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_y_continuous(expand = c(0, 0),
                                         breaks = seq_along(levels(df$Filename)),
                                         labels = levels(df$Filename)) +
                      theme_bw() +
                      theme(legend.position = "none",
                            panel.grid.minor = element_blank(),
                            panel.grid.major = element_blank()) +
                      labs(x = xLab, y = yLab)
                  
                  if (!is.null(userTheme)) scPlot <- scPlot + userTheme
                  
                  if (usePlotly){
                      scPlot <- scPlot +
                          theme(axis.ticks.y = element_blank(),
                                axis.text.y = element_blank(),
                                axis.title.y = element_blank(),
                                panel.grid = element_blank())
                      
                      status <- getSummary(x)
                      status <- status[status$Category == "Per base sequence content",]
                      status$Filename <- labels[status$Filename]
                      status$Filename <- factor(status$Filename, 
                                                levels = levels(df$Filename))
                      sideBar <- makeSidebar(status = status, key = key, 
                                             pwfCols = pwfCols)
                      
                      #plot dendro
                      if (dendrogram){
                          dx <- ggdendro::dendro_data(clusterDend)
                          dendro <- ggdend(dx$segments) +
                              coord_flip() +
                              scale_y_reverse(expand = c(0, 0)) +
                              scale_x_continuous(expand = c(0, 0.5))
                          dendro <- plotly::ggplotly(dendro, tooltip = NULL)
                          
                      }
                      else{
                          dendro <- plotly::plotly_empty()
                      }

                      
                      scPlot <- suppressWarnings(
                          suppressMessages(
                              plotly::subplot(dendro, sideBar, scPlot, 
                                              widths = c(0.1,0.08,0.82),
                                              margin = 0.001, shareY = TRUE) 
                          )
                      )
                      
                      ## This needs to be copied from master...
                      ## We have fixed it there...
                      ## manually edit tooltip to remove colour
                      sc <- lapply(seq_along(scPlot$x$data), function(x){
                          scPlot$x$data[[x]]$text <<- gsub("colour:.*<br />A", "A",  scPlot$x$data[[x]]$text)
                      })
                      
                  }
              }
              if (plotType == "line"){
                  df$Filename <- labels[df$Filename]
                  df <- df[!colnames(df) == "Base"]
                  df <- melt(df, id.vars = c("Filename", "Start"))
                  colnames(df) <- c("Filename", "Position", "Base", "Percent")
                  df$Base <- factor(df$Base, levels = c("T", "C", "A", "G"))
                  
                  #set colours
                  baseCols <- c(`T`="red", G = "black", A = "green", C = "blue")
                  
                  scPlot <- ggplot(df, aes_string(x = "Position", y = "Percent", colour = "Base")) +
                      geom_line() +
                      facet_wrap(~Filename) +
                      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                      scale_x_continuous(expand = c(0, 0)) +
                      guides(fill = FALSE) +
                      labs(x = "Position in read (bp)",
                           y = "Percent (%)") +
                      theme_bw() +
                      scale_colour_manual(values = baseCols)
                  
                  if (!is.null(userTheme)) scPlot <- scPlot + userTheme
                  
                  if(usePlotly){
                      scPlot <- suppressMessages(
                          ggplotly(scPlot) %>% layout(legend = list(x = 0.85, y = 1))
                      )
                  }
              }
              
              scPlot
          }
)
