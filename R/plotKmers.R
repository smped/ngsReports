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
#' # Not a great example
#' plotKmers(fdl[[1]])
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_bw theme
#' @importFrom ggplot2 xlab ylab labs
#' @importFrom magrittr %>%
#' @importFrom dplyr desc
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
            topK <- unique(df$Sequence[o])[1:n]
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
            refForX <- dplyr::data_frame(Base = as.character(refForX),
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
            if (n < length(pal)) pal <- pal[1:n]
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
            if (binned) kMerPlot <- kMerPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
            if (!is.null(userTheme)) kMerPlot <- kMerPlot + userTheme

            if(usePlotly){
              kMerPlot <- kMerPlot + ylab("Log2 Obs/Exp")
              kMerPlot <- suppressMessages(
                plotly::ggplotly(kMerPlot)
              )
            }
            kMerPlot
          }
)





# plotKmers <- function(x, subset, nKmers = 12, method = "overall",
#                             pwfCols, labels, naCol = "grey80",
#                             usePlotly = FALSE, clusterNames = FALSE,
#                             dendrogram = FALSE, plotType = "heatmap", nc = 2, ...){
#
#   # A basic cautionary check
#   stopifnot(grepl("(Fastqc|character)", class(x)))
#   stopifnot(is.numeric(nKmers))
#   stopifnot(method %in% c("overall", "individual"))
#
#   # Sort out the colours
#   if (missing(pwfCols)) pwfCols <- ngsReports::pwf
#   stopifnot(isValidPwf(pwfCols))
#   col <- getColours(pwfCols)
#
#
#   df <- tryCatch(Kmer_Content(x))
#
#   # Drop the suffix, or check the alternate labels
#   if (missing(labels)){
#     labels <- structure(gsub(".(fastq|fq|bam).*", "", fileName(x)),
#                         names = fileName(x))
#   }
#   else{
#     if (!all(fileName(x) %in% names(labels))) stop("All file names must be included as names in the vector of labels")
#   }
#   if (length(unique(labels)) != length(labels)) stop("The labels vector cannot contain repeated values")
#
#   # Get any arguments for dotArgs that have been set manually
#   dotArgs <- list(...)
#   if ("size" %in% names(dotArgs)){
#     lineWidth <- dotArgs$size
#   }
#   else{
#     lineWidth <- 0.2
#   }
#   if ("colour" %in% names(dotArgs) || "color" %in% names(dotArgs)){
#     i <- which(names(dotArgs) %in% c("colour", "color"))
#     lineCol <- dotArgs[[i]]
#   }
#   else{
#     lineCol <- "grey20"
#   }
#   allowed <- names(formals(ggplot2::theme))
#   keepArgs <- which(names(dotArgs) %in% allowed)
#   userTheme <- c()
#   if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])
#
#
#   # Get the top kMers
#   nKmers <- as.integer(nKmers)
#   if (method == "individual"){
#     topKmers <- split(df, f = df$Filename)
#     topKmers <- lapply(topKmers, function(x){
#       x <- dplyr::arrange(x, desc(x$Value))}) %>%
#       lapply(dplyr::slice, 1:nKmers) %>%
#       dplyr::bind_rows() %>%
#       magrittr::extract2("Sequence") %>%
#       unique()
#   }
#   if (method == "overall"){
#     topKmers <- dplyr::arrange(df, desc(df$Value)) %>%
#       dplyr::distinct(Sequence) %>%
#       dplyr::slice(1:nKmers) %>%
#       magrittr::extract2("Sequence")
#   }
#
#   df <- df[df$Sequence %in% topKmers,]
#
#   if(plotType == "heatmap"){
#
#   # Set the Sequence as a factor based on the first position it appears
#   # This way the colours will appear in order in the guide as well as the plot
#   kMerLevels <- dplyr::arrange(df, df$`Max_Obs/Exp_Position`) %>%
#     dplyr::distinct(Sequence) %>%
#     magrittr::extract2("Sequence")
#   df$Sequence <- factor(df$Sequence, levels = kMerLevels)
#
#   if (length(kMerLevels) < nKmers) {
#     message(paste("There is only data in the FASTQC reports for the top",
#                   length(kMerLevels),"Kmers."))
#     nKmers <- length(kMerLevels)
#   }
#
#   # Set the p-values to the -log10 scale
#   df$PValue <- -log10(df$PValue)
#   allInf <- all(is.infinite(df$PValue))
#   if (allInf) {
#     message(paste("All PValues for the top", nKmers, "Kmers are zero.\n",
#                   "This plot will be relatively uninformative."))
#     #if(usePlotly) stop(All kmer P-values equal zero)
#     df$PValue <- 0
#   }
#   else{
#     # If there are some finite -log10 transformed p-values,
#     # just add 1 to the the infinite ones
#     df$PValue[is.infinite(df$PValue)] <- max(df$PValue[is.finite(df$PValue)]) + 1
#   }
#
#
#   df <- reshape2::dcast(df, Filename~Sequence, value.var = "PValue")
#
#   if(clusterNames){
#     xx <- dplyr::select(df, -Filename)
#     xx[is.na(xx)] <- 0
#     clus <- as.dendrogram(hclust(dist(xx), method = "ward.D2"))
#     row.ord <- order.dendrogram(clus)
#     df <- df[row.ord,]
#   }
#   key <- df$Filename
#   df <- reshape2::melt(df, id.vars = "Filename", variable.name = "Sequence", value.name = "PValue")
#   df$Filename <- labels[df$Filename]
#   df$Filename <- factor(df$Filename, levels = unique(df$Filename))
#
#
#    # The inclusion of files without values needs to be worked on
#   # Maybe columns should be added during a dcast somewhere
#   if (allInf) {
#     # First change the missing values to ">0" as this is all the information we have
#     # This is done by gaining explicit NA values, then transforming
#
#     df <- dplyr::mutate(df, PValue = dplyr::if_else(is.na(PValue), ">0", "=0"))
#
#     kMerPlot <-  ggplot(df, aes_string(x = "Sequence", y = "Filename", fill = "PValue")) +
#       geom_tile(colour = "grey30", alpha = 0.9) +
#       scale_fill_manual(values = c(`=0` = getColours(pwfCols)["WARN"][[1]], `>0` = naCol)) +
#       labs(fill = "PValue")
#   }
#   else{
#     # First get explicit NA values
#
#     kMerPlot <- ggplot(df, aes_string(x = "Sequence", y = "Filename", fill = "PValue")) +
#       geom_tile(colour = "grey30", alpha = 0.9) +
#       scale_fill_gradient2(low = getColours(pwfCols)["PASS"],
#                            mid = getColours(pwfCols)["WARN"],
#                            high = getColours(pwfCols)["FAIL"],
#                            na.value = naCol,
#                            midpoint = max(df$PValue, na.rm = TRUE)/2) +
#       labs(fill = expression(paste(-log[10], "P")))
#   }
#
#   kMerPlot <- kMerPlot +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_y_discrete(expand = c(0, 0)) +
#     theme(panel.grid = element_blank(),
#           panel.background = element_blank())
#
#
#   if(usePlotly){
#     kMerPlot <- kMerPlot +
#       theme(axis.text = element_blank(), axis.ticks = element_blank(),
#             axis.title.y = element_blank()) +
#       labs(fill = "-log(10) P")
#
#     t <- getSummary(x)
#     t <- t[t$Category == "Kmer Content",]
#     t$Filename <- labels[t$Filename]
#     t$Filename <- factor(t$Filename, levels = unique(df$Filename))
#     t <- dplyr::right_join(t, unique(df["Filename"]), by = "Filename")
#
#
#     sideBar <- makeSidebar(t, pwfCols = pwfCols, key = key)
#
#
#     #plot dendrogram
#     if(dendrogram && clusterNames){
#
#       dx <- ggdendro::dendro_data(clus)
#       dendro <- ggdend(dx$segments) +
#         coord_flip() +
#         scale_y_reverse(expand = c(0, 0)) +
#         scale_x_continuous(expand = c(0,0.5))
#
#       kMerPlot <- plotly::subplot(dendro, sideBar, kMerPlot,
#                                   widths = c(0.1,0.08,0.82), margin = 0.001,
#                                   shareY = TRUE) %>%
#         plotly::layout(xaxis3 = list(title = "Kmer Sequence",
#                                      plot_bgcolor = "white"))
#     }else{
#
#       kMerPlot <- plotly::subplot(plotly::plotly_empty(),
#                                   sideBar,
#                                   kMerPlot,
#                                   widths = c(0.1,0.08,0.82),
#                                   margin = 0.001,
#                                   shareY = TRUE) %>%
#         plotly::layout(xaxis3 = list(title = "Kmer Sequence"),
#                        annotations = list(text = "Filename", showarrow = FALSE,
#                                           textangle = -90))
#   }
#
#   }else{
#     kMerPlot <- kMerPlot +
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#
#   }}
#   if(plotType == "line"){
#
#     df <- df %>%
#       dplyr::rename(Base = `Max_Obs/Exp_Position`) %>%
#       dplyr::mutate(Position = gsub("([0-9]*)-[0-9]*", "\\1", Base),
#                     Position = as.integer(Position)) %>%
#       dplyr::select(Filename, Sequence, Value, Position)   # The PValue and Count columns are ignored for this plot
#
#     # Set the x-axis for plotting
#     # As there will be significant gaps in this data,
#     # the bins used need to be obtained from another slot.
#     # The most complete will be Per_base_sequence_quality
#     # These values can then be incorporated in the final df for accurate plotting & labelling
#     refForX <- unique(Per_base_sequence_quality(x)$Base)
#     refForX <- dplyr::data_frame(Base = as.character(refForX),
#                                  Position = gsub("([0-9]*)-[0-9]*", "\\1", Base))
#     refForX$Position <- as.integer(refForX$Position)
#
#
#     # In order to get a line plot, zero values need to be added to the missing positions
#     # The above reference scale for X will be used to label the missing values
#     # Include all files to ensure all appear in the final plot
#     zeros <- with(df,
#                   expand.grid(list(Filename = Filename,
#                                    Sequence = unique(Sequence),
#                                    Value = 0,
#                                    Position = seq(0, max(Position) + 0.5, by = 0.5)),
#                               stringsAsFactors = FALSE))
#
#     # After the bind_rows, duplicate values will exist at some positions
#     # Spuriously introduced zeros need to be removed
#     df <- dplyr::bind_rows(df, zeros) %>%
#       dplyr::arrange(Filename, Sequence, Position, desc(Value)) %>%
#       dplyr::distinct(Filename, Sequence, Position, .keep_all = TRUE)
#
#     # Set the Sequence as a factor based on the first position it appears
#     # This way the colours will appear in order in the guide as well as the plot
#     kMerLevels <- df %>%
#       dplyr::filter(Value != 0) %>%
#       dplyr::arrange(Position) %>%
#       dplyr::distinct(Sequence) %>%
#       magrittr::extract2("Sequence")
#     df$Sequence <- factor(df$Sequence, levels = kMerLevels)
#
#     # Now draw the basic plots
#     kMerPlot <- ggplot(df,
#                        aes_string(x =" Position", y = "Value", colour = "Sequence")) +
#       geom_line() +
#       facet_wrap(~Filename, ncol = nc) +
#       scale_x_continuous(breaks = refForX$Position,
#                          labels = refForX$Base) +
#       theme_bw() +
#       ylab(expression(paste(log[2], " Obs/Exp"))) +
#       xlab("Position in read (bp)")
#
#     # Check for binned x-axis values to decied whether to rotate x-axis labels
#     # This should be clear if there are more than 2 characters in the plotted labels
#     binned <- max(nchar(dplyr::filter(refForX, Position %in% df$Position)$Base)) > 2
#     if (binned) {
#       kMerPlot <- kMerPlot +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#     }
#
#     if(usePlotly){
#       kMerPlot <- ggplotly(kMerPlot)
#     }
#   }
#   kMerPlot
# }
#
