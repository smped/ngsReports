#' @title Make the dedrogram for heatmap-style plots
#' 
#' @description Set the clusters for heatmap-style plots
#' 
#' @param df The data frame to be clustered
#' @param rowVal The rows to be clustered
#' @param colVal The value which will become column names
#' @param value The value to use for the clustering
#' 
#' @return A dendrogram
#' 
#' @examples 
#' # Get the files included with the package
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc", full.names = TRUE)[1:4]
#' cols <- c("Filename", "Position", "Illumina_Universal_Adapter")
#' ac <- Adapter_Content(fileList)[cols]
#' ngsReports:::makeDendrogram(ac, "Filename", "Position", "Illumina_Universal_Adapter")
#' 
#' @importFrom stats as.formula
#' 
#' @keywords internal
makeDendrogram <- function(df, rowVal, colVal, value){
    
    stopifnot(setequal(c(rowVal, colVal, value), names(df)))
    fm <- as.formula(paste0(rowVal, "~", colVal))
    mat <- reshape2::acast(df, fm, value.var = value)
    mat[is.na(mat)] <- 0
    clust <- hclust(dist(mat), method = "ward.D2")
    as.dendrogram(clust)
    
}
