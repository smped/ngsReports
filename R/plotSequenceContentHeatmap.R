#' @title Plot the per base content as a heatmap
#'
#' @description Plot the Per Base content for a set of FASTQC files.
#' Informative plot where per base sequence content (%A, %T, %G, %C),
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param subset \code{logical}. Return the values for a subset of files.
#' May be useful to only return totals from R1 files, or any other subset
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
#' @importFrom grDevices rgb
#'
#' @export

plotSequenceContent <- function(x, subset){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }

  # Get the NContent
  x <- tryCatch(x[subset])
  df <- tryCatch(Per_base_sequence_content(x))
  df <- Per_base_sequence_content(x)
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Base),
                      Start = as.integer(Start))
  df <- mutate(df, colour = rgb(rescale(A+0.33*G), rescale(T+0.33*G), rescale(C+0.33*G)))

  basicStat <- Basic_Statistics(x) %>% dplyr::select(Filename, Longest_sequence)

  df <- df %>% dplyr::right_join(basicStat, by = "Filename") %>%
    dplyr::select(Filename, Start, colour, Longest_sequence)

  dfInner <- df %>%
    split(f = .['Filename']) %>%
    lapply(function(x){
      dfFill <- data.frame(Start = 1:x$Longest_sequence[1])
      x <- dplyr::right_join(x, dfFill, by = "Start") %>%
        zoo::na.locf()
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(Start = as.integer(Start)) %>%
    dplyr::select(-Longest_sequence)

  sequenceContentHeatmap <- ggplot(dfInner,
                                            aes(x = Start,
                                                         y = Filename,
                                                         fill = colour)) +
    geom_tile() +
    scale_fill_manual(values = dfInner$colour) +
    theme(legend.position = "none",
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
  sequenceContentHeatmap
}

