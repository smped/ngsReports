#' @title Plot a summary of alignments
#'
#' @description Plot a summary of alignments from a set of log files
#'
#' @details Loads a set of alignment log files and creates a default plot.
#' Implemented aligners are \code{bowtie}, \code{bowtie2}, \code{Hisat2} and
#'  \code{STAR}.
#'
#' @param x Paths to one or more alignment log files
#' @param type The aligner used. Can be one of star, bowtie, bowtie2 or hisat2
#' @param usePlotly logical. If TRUE an interactive plot will be generated.
#' If FALSE a ggplot object will be output
#' @param ... Used to pass additional attributes to theme() and between methods
#' @param fill Colours used to fill the bars. Passed to scale_fill_manual.
#'
#' @return
#' A ggplot2 object, or a plotly object
#'
#' @examples
#'
#' f <- c("bowtie2PE.txt", "bowtie2SE.txt")
#' bowtie2Logs <- system.file("extdata", f, package = "ngsReports")
#' plotAlignmentSummary(bowtie2Logs, "bowtie2")
#'
#' @importFrom scales comma percent
#' @importFrom tidyselect one_of
#' @import ggplot2
#'
#' @export
plotAlignmentSummary <- function(
    x,
    type = c("star", "bowtie", "bowtie2", "hisat2"),
    usePlotly = FALSE,
    ...,
    fill = c("red", "yellow", "blue", rgb(0, 0.5, 1))
){
    ## Set the main arguments
    type <- match.arg(type)

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    ## Import the data
    df <- tryCatch(importNgsLogs(x, type))
    pFun <- paste0(".plot", stringr::str_to_title(type), "Alignment")
    args <- list(
        df = df,
        fill = fill
    )

    ## Generate the plot and add any theme information
    p <- do.call(pFun, args)
    if (!is.null(userTheme)) p <- p + userTheme

    if (usePlotly) {
        ## Set the tooltip
        tt <- c("x", "fill", "Percent", "Reads")
        ## Make interactive
        p <- plotly::ggplotly(
            p + theme(
                axis.title.y = element_blank(),
                legend.position = "none"
            ),
            tooltip = tt
        )
    }

    ## Return the plot
    p
}

.plotStarAlignment <- function(df, fill){

    cols <- c(
        "Filename",
        "Number_Of_Input_Reads",
        "Uniquely_Mapped_Reads_Number",
        "Number_Of_Reads_Mapped_To_Multiple_Loci",
        "Number_Of_Reads_Mapped_To_Too_Many_Loci"
    )
    df <- df[cols]
    subCols <- cols[seq(3, 5)]
    df$Unmapped_Other <- df$Number_Of_Input_Reads - rowSums(df[subCols])

    ## Remove the Log.final.out suffix
    df$Filename <- stringr::str_remove_all(df$Filename, "Log.final.out")
    ## Now gather for plotting
    df <- tidyr::gather(df, "Type", "Total", one_of(subCols))
    ## Remove the Percent/Number text from the Type column
    df$Type <-
        stringr::str_remove_all(df$Type, "(.+Of_Reads_|_Number$|_Percent$)")
    lv <-  c(
        "Uniquely_Mapped_Reads",
        "Mapped_To_Multiple_Loci",
        "Mapped_To_Too_Many_Loci",
        "Unmapped_Other"
    )
    df$Type <- factor(df$Type, levels = rev(lv))
    df$Percent = percent(df$Total / df$Number_Of_Input_Reads)
    df$Reads <- comma(df$Total)

    ggplot(
        df,
        aes_string(
            x = "Filename",
            y = "Total",
            fill = "Type",
            label1 = "Percent",
            label2 = "Reads"
        )
    ) +
        geom_bar(stat = "identity") +
        labs(y = "Total Reads") +
        scale_y_continuous(labels = comma, expand = expansion(c(0, 0.05))) +
        scale_fill_manual(values = fill) +
        coord_flip() +
        theme_bw()

}

.plotHisat2Alignment <- function(df, fill){


    df <- dplyr::select(df, -one_of("Paired_Reads", "Alignment_Rate"))
    keepCols <- names(df)[-seq_len(2)]
    df <- tidyr::gather(df, "Type", "Total", one_of(keepCols))
    df$Percent <- percent(df$Total / df$Total_Reads)
    df$Filename <- stringr::str_remove_all(df$Filename, ".(log|info|txt)$")
    df$Reads <- comma(df$Total)
    lv <-  c(
        "Unique_In_Pairs",
        "Unique_Unpaired",
        "Unique_Discordant_Pairs",
        "Multiple_Unpaired",
        "Multiple_In_Pairs",
        "Not_Aligned"
    )
    lv <- intersect(lv, df$Type)
    df$Type <- factor(df$Type, levels = rev(lv))
    fill <- grDevices::colorRampPalette(fill)(length(lv))

    ggplot(
        df,
        aes_string(
            x = "Filename",
            y = "Total",
            fill = "Type",
            label1 = "Percent",
            label2 = "Reads"
        )
    ) +
        geom_bar(stat = "identity") +
        labs(y = "Total Reads") +
        scale_y_continuous(labels = comma, expand = expansion(c(0, 0.05))) +
        scale_fill_manual(values = fill) +
        coord_flip() +
        theme_bw()

}

.plotBowtie2Alignment <- .plotHisat2Alignment

#' @importFrom tidyselect contains
.plotBowtieAlignment <- function(df, fill){

    df <- dplyr::select(df, -contains("Time"), -contains("Index"))
    keepCols <- names(df)[-seq_len(2)]
    df <- tidyr::gather(df, "Type", "Total", one_of(keepCols))
    df <- dplyr::filter(df, !is.na(Total))
    df$Percent <- percent(df$Total / df$Reads_Processed)
    df$Filename <- stringr::str_remove_all(df$Filename, ".(log|info|txt)$")
    df$Reads <- comma(df$Total)

    df$Type <- factor(df$Type, levels = rev(unique(df$Type)))
    nLevels <- length(levels(df$Type))
    if (nLevels <= length(fill)) {
        fill <- fill[seq_len(nLevels)]
    }
    else {
        fill <- grDevices::colorRampPalette(fill)(nLevels)
    }

    ggplot(
        df,
        aes_string(
            x = "Filename",
            y = "Total",
            fill = "Type",
            label1 = "Percent",
            label2 = "Reads"
        )
    ) +
        geom_bar(stat = "identity") +
        labs(y = "Total Reads") +
        scale_y_continuous(labels = comma, expand = expansion(c(0, 0.05))) +
        scale_fill_manual(values = fill) +
        coord_flip() +
        theme_bw()

}