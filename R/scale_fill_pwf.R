# #' @title Construct a scale using PwfCols
# #'
# #' @description Construct a scale using PwfCols
# #'
# #' @details This constructs a ggplot scale using the values contained in the fill aesthetic
# #' and the supplied breaks for PASS/WARN/FAIL criteria.
# #' As this doesn't follow the conventional ggplot syntax using more of functional approach,
# #' it will be a hidden function
# #' @include PwfCols.R
# #' @importFrom ggplot2 scale_fill_gradientn
scale_fill_pwf <- function(vals, pwfCols, breaks = c(0, 5, 10, 100), passLow = TRUE){

  gradCols <- getColours(pwfCols) # Get the default colours
  upr <- suppressWarnings(max(vals, na.rm = TRUE)) # The maximum colour in the data
  nCols <- dplyr::if_else(upr == 0, 1, findInterval(upr, breaks) + 1 ) # The number of breaks
  # Get the gradient based on the upper & lower breakpoints around the max value
  # Then define the upper colour based on where `upr` lies in this range
  upperGradient <- grDevices::colorRampPalette(c(gradCols[max(1, nCols - 1)], gradCols[nCols]))(100)
  maxCol <- findInterval(upr, seq(breaks[max(1, nCols - 1)], breaks[nCols], length.out = 100))
  maxCol <- upperGradient[maxCol]
  nCols <- dplyr::if_else(maxCol == gradCols[1], 1, nCols)
  if (nCols == 1){
    ggplot2::scale_fill_gradient(low = gradCols[1], high = gradCols[1], na.value = gradCols[1])
  }
  else{
    gradCols[nCols] <- maxCol
    breaks[nCols] <- upr
    if(passLow) ggplot2::scale_fill_gradientn(colours = gradCols[1:3], values = breaks / upr, na.value = gradCols[1])
    else ggplot2::scale_fill_gradientn(colours = gradCols[3:1], values = breaks / upr, na.value = gradCols[1])

    }
}
