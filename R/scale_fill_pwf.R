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
  
  # # Find the number of categories present in the data
  # rng <- suppressWarnings(range(vals, na.rm = TRUE))
  # minCat <- findInterval(rng[1], breaks)
  # maxCat <- findInterval(rng[2], breaks)
  # 
  # # If only 1 category is present, find which category & create a gradient
  # if (minCat == maxCat){
  #   # Create a gradient & find the colour limits
  #   fullGradient <- grDevices::colorRampPalette(c(gradCols[minCat], gradCols[minCat + 1]))(100)
  #   indexVec <- seq(breaks[minCat], breaks[minCat + 1])
  #   minCol <- fullGradient[findInterval(rng[1], indexVec)]
  #   maxCol <- fullGradient[findInterval(rng[2], indexVec)]
  #   if(missing(na.value)) na.value <- gradCols[["MAX"]]
  #   return(
  #     ggplot2::scale_fill_gradient(low = minCol, high = maxCol, na.value = na.value)
  #   )
  # }
  # 
  # # If there is more than one category:
  # # The upper and lower gradients need to be set separately to find the appropriate extrema
  # lowerGradient <- grDevices::colorRampPalette(c(gradCols[minCat], gradCols[minCat + 1]))(100)
  # indexVec <- seq(breaks[minCat], breaks[minCat + 1], length.out = 100)
  # minCol <- lowerGradient[findInterval(rng[1], indexVec)]
  # upperGradient <- grDevices::colorRampPalette(c(gradCols[maxCat], gradCols[maxCat + 1]))(100)
  # indexVec <- seq(breaks[maxCat], breaks[maxCat + 1], length.out = 100)
  # maxCol <- upperGradient[findInterval(rng[2], indexVec)]
  # # Now reset the colours to reflect the new extrema
  # gradCols[minCat] <- minCol
  # gradCols[maxCat + 1] <- maxCol
  # # Reset the breaks
  # breaks[minCat] <- rng[1]
  # breaks[maxCat + 1] <- rng[2]
  # # Reset the breaks & basic colours to match the number of categories
  # gradCols <- gradCols[seq(minCat, maxCat + 1)]
  # breaks <- breaks[seq(minCat, maxCat + 1)]
  # 
  # if(missing(na.value)) na.value <- gradCols[1]
  # ggplot2::scale_fill_gradientn(colours = gradCols, values = scales::rescale(breaks),
  #                               limits = rng, na.value = na.value)

  # # We need to set the upper gradient & the lower gradient
  #
  # # upr <- suppressWarnings(max(vals, na.rm = TRUE)) # The maximum value in the data
  # nCols <- dplyr::if_else(upr == 0, 1, findInterval(upr, breaks) + 1 ) # The number of breaks
  # # Get the gradient based on the upper & lower breakpoints around the max value
  # # Then define the upper colour based on where `upr` lies in this range
  # upperGradient <- grDevices::colorRampPalette(c(gradCols[max(1, nCols - 1)], gradCols[nCols]))(100)
  # maxCol <- findInterval(upr, seq(breaks[max(1, nCols - 1)], breaks[nCols], length.out = 100))
  # maxCol <- upperGradient[maxCol]
  # nCols <- dplyr::if_else(maxCol == gradCols[1], 1, nCols)
  # if (nCols == 1){
  #   ggplot2::scale_fill_gradient(low = gradCols[1], high = gradCols[1], na.value = gradCols[1])
  # }
  # else{
  #   gradCols[nCols] <- maxCol
  #   breaks[nCols] <- upr
  #   if(passLow) {
  #     o <- seq(1, nCols, by = 1)
  #   }
  #   else{
  #     o <- seq(nCols, 1, by = -1)
  #   }
  #   ggplot2::scale_fill_gradientn(colours = gradCols[o], values = breaks[1:nCols] / upr,
  #                                 na.value = gradCols[1])
  #   }
}
