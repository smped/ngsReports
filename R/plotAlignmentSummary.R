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
#'
#' @return
#' A ggplot2 object, or a plotly object
#'
#' @importFrom scales comma percent
#' @import ggplot2
#'
#' @export
plotAlignmentSummary <- function(
    x,
    type = c("star", "bowtie", "bowtie2", "hisat2"),
    usePlotly = FALSE
){
    ## Set the main arguments
    type <- match.arg(type)

    ## Set the tooltip
    if (type == "star") {
        tt <- c("x", "fill", "Percent", "Reads")
    }
    else{
        stop("not done yet")
    }

    ## Import the data
    df <- tryCatch(importNgsLogs(x, type))
    pFun <- paste0(".plot", str_to_title(type), "Alignment")
    args <- list(
        df = df
    )
    p <- do.call(pFun, args)
    if (usePlotly) {
        p <- plotly::ggplotly(
            p + theme(legend.position = "none"),
            tooltip = tt
        )
    }

    ## And return the plot
    p
}

.plotStarAlignment <- function(
    df,
    fill =rgb(c(1, 1, 0, 0), c(0, 1, 0, 0.5), c(0, 0, 1, 1))
){

    cols <- c(
        "Filename",
        "Number_Of_Input_Reads",
        "Uniquely_Mapped_Reads_Number",
        "Number_Of_Reads_Mapped_To_Multiple_Loci",
        "Number_Of_Reads_Mapped_To_Too_Many_Loci"
    )
    df <- df[cols]
    df$Unmapped_Other <- df$Number_Of_Input_Reads - rowSums(df[,3:5])

    ## Remove the Log.final.out suffix
    df$Filename <- stringr::str_remove_all(df$Filename, ".Log.final.out")
    ## Now gather for plotting
    df <- tidyr::gather(df, "Type", "Total", -(1:2))
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
        scale_y_continuous(labels = comma, expand = expand_scale(c(0, 0.05))) +
        scale_fill_manual(values = fill) +
        coord_flip() +
        theme_bw()

}
