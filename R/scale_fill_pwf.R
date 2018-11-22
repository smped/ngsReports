#' @title Construct a scale using PwfCols
#'
#' @description Construct a scale using PwfCols
#'
#' @details This constructs a ggplot scale using the values contained in the fill aesthetic
#' and the supplied breaks for PASS/WARN/FAIL criteria.
#' As this doesn't follow the conventional ggplot syntax, using more of functional approach,
#' it will be a hidden function
#' 
#' @param vals The values which need to have the scale generated
#' @param pwfCols An object of class PwfCols
#' @param breaks The breaks for the PWF bins
#' @param passLow Is the PASS category at the low or high end of the numeric range
#' @param na.value The colour to plot for missing values
#' 
#' @return 
#' Returns a ggplot scale object
#'
#' @include PwfCols.R
#' @import ggplot2
#'
#' @keywords internal
#'
scale_fill_pwf <- function(vals, pwfCols, breaks = c(0, 5, 10, 100), passLow = TRUE, na.value = "white"){
    
    # passLow defines whether pass is the low score or the high score
    # organise the colours based on this
    o <- seq(1, 4)
    if (!passLow) o <- rev(o)
    gradCols <- getColours(pwfCols)[o] # Get the default gradient colours
    
    # Find which of the pwf bins are present
    bins <- cut(vals, breaks = breaks, include.lowest = TRUE)
    bins <- range(as.integer(bins))
    bins <- unique(bins)
    # Create an even sequence between the min & max of the range
    n <- seq(breaks[min(bins)], breaks[max(bins) + 1], length.out = 101)
    n <- cut(n, breaks = breaks)
    n <- split(n, n)
    # Now create a colour vector between the extreme points
    cols <- lapply(bins, function(x){
        l <- length(n[[x]]) + 1
        colorRampPalette(gradCols[c(x,x + 1)])(l)[-l]
    })
    cols <- as.character(c(unlist(cols), gradCols[max(bins) + 1]))
    # Remove any breaks outside of the range
    breaks <- breaks[seq(min(bins), max(bins) + 1)]
    # Return the scale object
    scale_fill_gradientn(colours = cols, breaks = breaks, limits = range(breaks), na.value = na.value)
    
}
