#' @title Draw a PCA plot for Fast QC modules
#'
#' @description Draw a PCA plot for Fast QC modules across multiple samples
#'
#' @details
#' This carries out PCA on all or a subset of FastQC modules and plots the
#' output using either ggplot or plotly. Clustering of the PCA can be carried
#' out using a k-means approach.
#'
#'
#' @param x Can be a \code{FastqcData}, \code{FastqcDataList} or file paths
#' @param module \code{character} vector containing 
#'  the desired FastQC module (eg. c("Per_base_sequence_quality","Per_base_sequence_content"))
#' @param usePlotly \code{logical}. Output as ggplot2 (default) or plotly
#' object.
#' @param pwfCols Object of class \code{\link{PwfCols}} containing the colours
#' for PASS/WARN/FAIL
#' @param labels An optional named vector of labels for the file names.
#' All filenames must be present in the names.
#' File extensions are dropped by default
#' @param cluster \code{logical} default \code{FALSE}. If set to \code{TRUE},
#' fastqc data will be clustered using hierarchical clustering
#' @param ... Used to pass additional attributes to theme() and between methods
#'
#' @return A standard ggplot2 object, or an interactive plotly object
#'
#' @examples
#'
#'
#' @docType methods
#'
#' @importFrom dplyr left_join
#' @import ggplot2
#' @import mclust 
#' @import tibble
#' 
#'
#' @name plotFastqcPCA
#' @rdname plotFastqcPCA-methods
#' @export
setGeneric("plotFastqcPCA", function(
    x, usePlotly = FALSE, labels, pwfCols, plotValue, cluster = FALSE, ...){
    standardGeneric("plotFastqcPCA")
}
)
#' @rdname plotFastqcPCA-methods
#' @export
setMethod("plotFastqcPCA", signature = "ANY", function(
    x, usePlotly = FALSE, labels, pwfCols, plotValue, cluster = FALSE, ...){
    .errNotImp(x)
}
)
#' @rdname plotFastqcPCA-methods
#' @export
setMethod("plotFastqcPCA", signature = "character", function(
    x, usePlotly = FALSE, labels, pwfCols, plotValue, cluster = FALSE, ...){
    x <- FastqcDataList(x)
    if (length(x) == 1) x <- x[[1]]
    plotFastqcPCA(x, usePlotly, labels, plotValue, pwfCols, cluster, ...)
}
)
#' @rdname plotFastqcPCA-methods
#' @export
setMethod("plotFastqcPCA", signature = "FastqcDataList", function(
    x, usePlotly = FALSE, labels, pwfCols, plotValue, module, cluster = FALSE, ...){
    
    # if(modules == "all") modules <- c("Per_base_sequence_quality", "Per_tile_sequence_quality", 
    #                                    "Per_sequence_quality_scores", "Per_base_sequence_content", 
    #                                    "Per_sequence_GC_content", "Per_base_N_content", 
    #                                    "Sequence_Length_Distribution", "Sequence_Duplication_Levels", 
    #                                    "Overrepresented_sequences", "Adapter_Content", "Kmer_Content", 
    #                                    "Total_Deduplicated_Percentage")
    
    
    pFun <- paste0(".generate", str_to_title(module), "PCA")
    args <- list(
        x = x)
    
    df <- do.call(pFun, args)
    pca <- prcomp(df, scale. = TRUE)
    
    scores <- as_tibble(pca$x) 
    
    dis <- scores[1:2]
    Distance <- dist(dis)
    
    
    #https://www.statmethods.net/advstats/cluster.html
    k <- Mclust(scores[1:2])$G
    kM <- kmeans(Distance, centers = k)
    
    
    kM <- kM$cluster 
    
    splitClus <- split(names(kM), kM)
    
    clusterDF <- lapply(1:length(splitClus), function(x){
        
        data.frame(Filename = splitClus[[x]], cluster = as.character(x))
        
        
    }) 
    
    clusterDF <- bind_rows(clusterDF)
    
    
    scores <- rownames_to_column(scores, "Filename")
    
    data <- left_join(scores, clusterDF, by = "Filename")
    
    
    ## get convex edges
    find_hull <- function(data) data[chull(data$PC2, data$PC1), ]
    
    hulls <- plyr::ddply(data, "cluster", find_hull)
    hulls$cluster <- factor(hulls$cluster, levels = unique(hulls$cluster))
    
    
    p <- ggplot(data = data) +
        geom_polygon(data = hulls, aes(x = PC1, y = PC2, fill = cluster), alpha = 0.2) +
        geom_point(aes(group = Filename, x = PC1, y = PC2)) +
        geom_hline(yintercept=0, colour="darkgrey") + 
        geom_vline(xintercept=0, colour="darkgrey") 
    
    plotlyPlot <- ggplotly(p)
    
    s <- split(data, data$cluster)
    
    plotlyPlot$x$data[1:k] <- lapply(1:k, function(j){
        
        names <- s[[j]]$Filename
        names <- paste(names, collapse = "<br>")
        
        plotlyPlot$x$data[[j]]$text <- names
        plotlyPlot$x$data[[j]]
    })
    
    plotlyPlot
}
)


.generatePer_base_sequence_qualityPCA <- function(x){
    
    df <- getModule(x, "Per_base_sequence_quality")
    df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))
    
    
    ## Adjust the data for files with varying read lengths
    ## This will fill NA values with the previous values
    df <- lapply(split(df, f = df$Filename), function(y){
        Longest_sequence <-
            gsub(".*-([0-9]*)", "\\1", as.character(y$Base))
        Longest_sequence <- max(as.integer(Longest_sequence))
        dfFill <- data.frame(Start = seq_len(Longest_sequence))
        y <- dplyr::right_join(y, dfFill, by = "Start")
        na.locf(y)
    })
    
    df <- dplyr::bind_rows(df)[c("Filename", "Start", "Mean")]
    
    df <- dcast(df, Filename ~ factor(as.character(df$Start), 
                                      levels = unique(as.character(df$Start))), 
                value.var = "Mean", fill = 0) 
    df <- column_to_rownames(df, "Filename")
    
    df 
}

#c("Per_base_sequence_quality", "Per_tile_sequence_quality", 
#                                    "Per_sequence_quality_scores", "Per_base_sequence_content", 
#                                    "Per_sequence_GC_content", "Per_base_N_content", 
#                                    "Sequence_Length_Distribution", "Sequence_Duplication_Levels", 
#                                    "Overrepresented_sequences", "Adapter_Content", "Kmer_Content", 
#                                    "Total_Deduplicated_Percentage")


.generatePer_sequence_quality_scoresPCA <- function(x){
    
    df <- getModule(x, "Per_sequence_quality_scores")
    
    df <- dcast(df, Filename ~ factor(as.character(df$Quality), 
                                      levels = unique(as.character(df$Quality))), 
                value.var = "Count", fill = 0) 
    df <- column_to_rownames(df, "Filename")
    
    df 
}  


.generatePer_sequence_GC_contentPCA <- function(x){
    
    df <- getModule(x, "Per_sequence_GC_content")
    
    df <- dcast(df, Filename ~ factor(as.character(df$GC_Content), 
                                      levels = unique(as.character(df$GC_Content))), 
                value.var = "Count", fill = 0) 
    df <- column_to_rownames(df, "Filename")
    
    df 
}  


.generatePer_base_sequence_contentPCA <- function(x){
    
    df <- getModule(x, "Per_base_sequence_content")
    
    df$Start <- as.integer(gsub("([0-9]*)-[0-9]*", "\\1", df$Base))
    
    
    ## Adjust the data for files with varying read lengths
    ## This will fill NA values with the previous values
    df <- lapply(split(df, f = df$Filename), function(y){
        Longest_sequence <-
            gsub(".*-([0-9]*)", "\\1", as.character(y$Base))
        Longest_sequence <- max(as.integer(Longest_sequence))
        dfFill <- data.frame(Start = seq_len(Longest_sequence))
        y <- dplyr::right_join(y, dfFill, by = "Start")
        na.locf(y)
    })
    
    df <- dplyr::bind_rows(df)[c("Filename", "Start", "G", "A", "T", "C")]
    
    t <- lapply(c("G", "A", "T", "C"), function(x){
        df <- dcast(df, Filename ~ factor(as.character(df$Start), 
                                          levels = unique(as.character(df$Start))), 
                    value.var = x, fill = 0) 
        
        colnames(df)[2:ncol(df)] <- paste0(colnames(df)[2:ncol(df)], ".", x)
        df
    })
    
    df <- Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="Filename"), t)
}  

.generateSequence_Length_DistributionPCA <- function(x){
    
    df <- getModule(x, "Sequence_Length_Distribution")
    
    df <- lapply(split(df, f = df$Filename), function(y){
        Longest_sequence <-
            gsub(".*-([0-9]*)", "\\1", as.character(y$Length))
        Longest_sequence <- max(as.integer(Longest_sequence))
        dfFill <- data.frame(Lower = seq_len(Longest_sequence))
        y <- dplyr::right_join(y, dfFill, by = "Lower")
        na.locf(y)
    })
    
    df <- dplyr::bind_rows(df)[c("Filename", "Lower", "Count")]
    
    df <- dcast(df, Filename ~ factor(as.character(df$Lower), 
                                      levels = unique(as.character(df$Lower))), 
                value.var = "Count", fill = 0) 
    df <- column_to_rownames(df, "Filename")
    
    df 
} 


