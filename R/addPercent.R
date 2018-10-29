#' @title Add a percentage sign
#' 
#' @description Add a percentage sign to the end of a string
#' 
#' @param x Any vector
#' 
#' @return character vector
#' 
#' @examples 
#' 
#' x <- 1:10
#' ngsReports:::addPercent(x)
#' 
#' @keywords internal
#' 
addPercent <- function(x){
    if (is.factor(x)) message("Factors will be converted to characters")
    paste0(x, "%")
}