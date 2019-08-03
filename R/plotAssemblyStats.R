#' @title Plot a summary of assembly logs
#'
#' @description Plot a summary of assembly stats from a set of log files
#'
#' @details Loads a set of assembly log files and creates a default plot.
#' Implemented tools are \code{quast} and \code{BUSCO}.
#' quast will plot a parralel coordinate plot of some assembly statistics
#' BUSCO will plot a stacked barplot of completeness statistics
#'
#' @param x Paths to one or more log files
#' @param type The tool used. Can be one of quast or busco
#' @param usePlotly logical. If TRUE an interactive plot will be generated.
#' If FALSE a ggplot object will be output
#' @param plotType \code{character}. Plot type to output, one of bar or paracoord.
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return
#' A ggplot2 object, or a plotly object
#'
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#'
#' @export
plotAssemblyStats <- function(x, type = "quast", usePlotly = FALSE, plotType = "bar", ...){
    
    ## Set the main arguments
    type <- match.arg(type)
    
    ## Get any arguments for dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(theme))
    keepArgs <- which(names(dotArgs) %in% allowed)
    userTheme <- c()
    if (length(keepArgs) > 0) userTheme <- do.call(theme, 
                                                   dotArgs[keepArgs])
    
    ## Import the data
    df <- tryCatch(importNgsLogs(x, 
                                 type))
    pFun <- paste0(".plot", 
                   str_to_title(type), "Stats")
    args <- list(
        df = df,
        usePlotly = usePlotly,
        plotType = plotType
    )
    
    ## Generate the plot and add any theme information
    p <- do.call(pFun, 
                 args)
    # if (!is.null(userTheme)) p <- p + userTheme
    
    ## Return the plot
    p
    
    
    
    
    
}

.plotQuastStats <- function(df, usePlotly, plotType){
    
    #df <- importNgsLogs(x, type = "quast")
    labels <- .makeLabels(df, 
                                       pattern = ".tsv", 
                                       col = "fileNames")
    df$fileNames <- labels[df$fileNames]
    
    
    dfLong <- tidyr::gather(df,
                            "variable", 
                            "Value", 
                            2:7)
    varaible <- c()
    
    dfLong$fileNames <- factor(dfLong$fileNames, 
                               levels = unique(dfLong$fileNames))
    dfLong$variable <- factor(dfLong$variable,
                              levels = unique(dfLong$variable))
    
    
    
    
    
    # data <- ggparcoord(df, columns = 2:7, groupColumn = 1, showPoints = TRUE, scale = "uniminmax")$data
    # 
    # data <- dplyr::left_join(data, df2, by = c("fileNames", "variable"))
    
    if(plotType == "parcoord"){
        
        
        dfLong <- dplyr::group_by(dfLong, 
                                  variable)
        dfLong <- dplyr::mutate(dfLong, 
                                scaled = (Value-min(Value))/(max(Value)-min(Value)))
        dfLong$Value <- paste(round(dfLong$Value/1000000, 2), "Mb")
        p <- suppressWarnings(ggplot() + 
                                  geom_line(data = dfLong,
                                            aes_string(x = "variable",
                                                       y = "scaled", 
                                                       colour = "fileNames", 
                                                       group = "fileNames"),
                                            size = 1, linetype = 2) +
                                  geom_point(data = dfLong, aes_string(x = "variable",
                                                                       y = "scaled",
                                                                       colour = "fileNames",
                                                                       group = "fileNames", 
                                                                       Value = "Value"), 
                                             size = 4, shape = "diamond") +
                                  theme(axis.title.y = element_blank(), 
                                        axis.ticks.y = element_blank(), 
                                        axis.text.y = element_blank(), 
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.major.x = element_line(colour = "black")) +
                                  scale_colour_brewer(palette = "Dark2") + 
                                  scale_y_continuous(expand = c(0.1,0)) + xlab("Metric"))
        if(usePlotly) p <- ggplotly(p, tooltip = c("Value", "group"))
        else{
            # set up min max values for axis annotations 
            minMaxL <- lapply(c("min", "max"), function(pass){
                sapply(colnames(df)[2:7],  function(y){
                    n <- do.call(pass, list(df[[y]])) 
                    if(n > 100000) paste(round(n/1000000, 2), "Mb")
                    else paste(round(n/1000, 2), "Kb")
                })
            }) 
            minDf <- tibble(x = rep(1:6 + 0.04, 2), 
                            y = rep(0:1, each = 6), 
                            label = unlist(minMaxL),
                            hjust = rep(-0.05, 12),
                            vjust = rep(c(1.5,-0.5), each = 6))
            
            p <- p + geom_text(data = minDf,
                               aes_string(x = "x", 
                                          y = "y", 
                                          label = "label", 
                                          vjust = "vjust",
                                          hjust = "hjust"))
        }
        
    } 
    else{
        
        p <- ggplot() + 
            geom_bar(data = dfLong,
                     aes_string(x = "fileNames",
                                y = "Value", 
                                fill = "fileNames",
                                group = "fileNames"),
                     stat = "identity")  +
            facet_wrap(~variable, 
                       scales = "free_y") + 
            theme_linedraw() +
            scale_fill_brewer(palette = "Dark2")
        
    }
    
    p
    
    
    
}

.plotBuscoStats <- function(df, usePlotly){
    
    #df <- importNgsLogs(x, type = "busco")
    
    df <- tidyr::gather(df, 
                        "Status", 
                        "Count", 
                        2:5)
    df[["name"]] <- factor(df[["name"]],
                           levels = unique(df[["name"]]))
    df[["Status"]] <- factor(df[["Status"]], 
                             levels = c("missing", 
                                        "fragmented", 
                                        "completeDuplicated",
                                        "completeSingleCopy"))
    
    p <- ggplot() + 
        geom_bar(data = df, 
                 aes_string(x = "name", 
                            y = "Count", 
                            fill = "Status", 
                            group = "Status"),
                 stat = "identity", width = 0.8) + 
        coord_flip() +
        scale_fill_manual(values = c(completeSingleCopy = "darkgreen", 
                                     completeDuplicated = "yellow", 
                                     fragmented = "darkorange", 
                                     missing = "darkred")) + 
        theme_linedraw()
    if(usePlotly) p <- ggplotly(p, 
                                tooltip = c("Count", "group"))
    
    p
}
