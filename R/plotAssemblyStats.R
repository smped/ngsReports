#' @title Plot a summary of assembly logs
#'
#' @description Plot a summary of assembly stats from a set of log files
#'
#' @details Loads a set of assembly log files and creates a default plot.
#' Implemented tools are `quast` and `BUSCO`.
#' quast will plot a parralel coordinate plot of some assembly statistics
#' BUSCO will plot a stacked barplot of completeness statistics
#'
#' @param x Paths to one or more log files
#' @param type The tool used. Can be one of quast or busco
#' @param usePlotly logical. If TRUE an interactive plot will be generated.
#' If FALSE a ggplot object will be output
#' @param plotType `character`. Plot type to output, one of bar or
#' paracoord.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return
#' A ggplot2 object, or a plotly object
#'
#' @examples
#'
#' #getquast log filenames
#' quastFiles <- system.file("extdata",
#' c("quast1.tsv", "quast2.tsv"), package = "ngsReports")
#'
#' # The default plot
#' plotAssemblyStats(quastFiles)
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom plotly ggplotly
#'
#' @export
plotAssemblyStats <- function(
        x, type = c("quast", "busco"), usePlotly = FALSE,
        plotType = c("bar", "paracoord"), ...){

    ## Set the main arguments
    type <- match.arg(type)
    plotType <- match.arg(plotType)

    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, dotArgs[keepArgs])

    ## Import the data
    df <- tryCatch(importNgsLogs(x,  type))
    pFun <- paste0(".plot", str_to_title(type), "Stats")
    args <- list(
        df = df,
        usePlotly = usePlotly,
        plotType = plotType
        ###############################################################
        ## dotArgs need to go here to handle themes within each plot ##
        ###############################################################
    )

    ## Generate the plot and add any theme information
    p <- do.call(pFun, args)

    ## Return the plot
    p

}

.plotQuastStats <- function(df, usePlotly, plotType, ...){

    labels <- .makeLabels(df, pattern = ".tsv", col = "fileNames")
    df$fileNames <- labels[df$fileNames]

    dfLong <- tidyr::gather(df, "variable", "Value", seq(2, 7))
    variable <- Value <- c()

    dfLong$fileNames <- factor(
        dfLong$fileNames,
        levels = unique(dfLong$fileNames)
    )
    dfLong$variable <- factor(
        dfLong$variable,
        levels = unique(dfLong$variable)
    )
    variable <- scaled <- fileNames <- c()

    if (plotType == "paracoord") {

        dfLong <- dplyr::group_by(dfLong, variable)
        dfLong <- dplyr::mutate(
            dfLong,
            scaled = (Value - min(Value)) / (max(Value) - min(Value))
        )
        dfLong$Value <- paste(round(dfLong$Value/1000000, 2), "Mb")
        p <- suppressWarnings(
            ggplot() +
                geom_line(
                    data = dfLong,
                    aes(
                        variable, scaled, colour = fileNames, group = fileNames
                    ),
                    linewidth = 1, linetype = 2
                ) +
                geom_point(
                    data = dfLong,
                    aes(
                        variable, scaled, colour = fileNames, group = fileNames,
                        Value = Value
                    ),
                    size = 4,
                    shape = "diamond"
                ) +
                theme(
                    axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank(),
                    legend.key = element_rect(fill = "white"),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.major.x = element_line(colour = "black")
                ) +
                scale_colour_brewer(palette = "Dark2") +
                scale_y_continuous(expand = c(0.1,0)) + xlab("Metric")
        )
        if (usePlotly) p <- ggplotly(p, tooltip = c("Value", "group"))
        else{
            # set up min max values for axis annotations
            minMaxL <- lapply(
                c("min", "max"),
                FUN = function(pass){
                    vapply(
                        colnames(df)[seq(2, 7, by = 1)],
                        FUN.VALUE = character(1),
                        FUN = function(y){
                            n <- do.call(pass, list(df[[y]]))
                            if(y %in% colnames(df)[seq(2, 5, by = 1)]){
                                if (n > 100000) paste(round(n/1000000, 2), "Mb")
                                else paste(round(n/1000, 2), "Kb")
                            }
                            else as.character(n)

                        }
                    )
                }
            )
            minDf <- tibble(
                x = rep(seq_len(6) + 0.04, 2),
                y = rep(c(0, 1), each = 6),
                label = unlist(minMaxL),
                hjust = rep(-0.05, 12),
                vjust = rep(c(1.5,-0.5), each = 6)
            )

            label <- vjust <- hjust <- c()
            p <- p +
                geom_text(
                    data = minDf,
                    aes(x, y, label = label, vjust = vjust, hjust = hjust)
                )
        }

    }
    else{

        p <- ggplot() +
            geom_bar(
                data = dfLong,
                aes(fileNames, Value, fill = fileNames, group = fileNames),
                stat = "identity"
            )  +
            facet_wrap(~variable, scales = "free_y") +
            theme_bw() +
            scale_fill_brewer(palette = "Dark2")

    }

    p

}

.plotBuscoStats <- function(df, usePlotly, ...){

    df <- tidyr::gather(df, "Status", "Count", seq(2, 5))
    df[["name"]] <- factor(df[["name"]], levels = unique(df[["name"]]))
    df[["Status"]] <- factor(
        df[["Status"]],
        levels = c(
            "missing", "fragmented", "completeDuplicated", "completeSingleCopy"
        )
    )

    name <- c()
    p <- ggplot() +
        geom_bar(
            data = df, aes(name, Count, fill = Status, group = Status),
            stat = "identity", width = 0.8
        ) +
        coord_flip() +
        scale_fill_manual(
            values = c(
                completeSingleCopy = "darkgreen",
                completeDuplicated = "yellow",
                fragmented = "darkorange",
                missing = "darkred"
            )
        ) +
        theme_bw()
    if (usePlotly) p <- ggplotly(p,  tooltip = c("Count", "group"))

    p
}
