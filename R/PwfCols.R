#' The PwfCols class and associated methods
#'
#' @description Define the PwfCols class and associated methods
#'
#' @details
#' This is an object of with four colours in components named PASS, WARN, FAIL and MAX.
#' Used to indicate these categories as defined on the standard plots from fastqc
#'
#' @slot PASS A vector of length 1, defining the colour for PASS in rgb format. Defaults to rgb(0, 0.8, 0)
#' @slot WARN A vector of length 1, defining the colour for WARN in rgb format. Defaults to rgb(0.9, 0.9, 0.2)
#' @slot FAIL A vector of length 1, defining the colour for FAIL in rgb format. Defaults to rgb(0.8, 0.2, 0.2)
#' @slot MAX A vector of length 1, defining the colour for an extreme FAIL or NA in rgb format.
#' Defaults to rgb(1, 1, 1)
#'
#' @include validationFunctions.R
#' @aliases PwfCols
setClass("PwfCols", slots = c(PASS = "character",
                              WARN = "character",
                              FAIL = "character",
                              MAX = "character"))
setValidity("PwfCols", isValidPwf)

#' @title Get and Set the Colours for a PwfCols Objects
#'
#' @description Define the generics & methods for setting/getting colours from an object of class PwfCols
#'
#' @param object An object of class PwfCols
#'
#' @include AllGenerics.R
#'
#' @export
#' @rdname getColours
#' @aliases getColours
setMethod("getColours", "PwfCols", function(object){
  vals <- c(object@PASS,
            object@WARN,
            object@FAIL,
            object@MAX)
  names(vals) <- c("PASS", "WARN", "FAIL", "MAX")
  vals
})

#' @export
#' @aliases getColors
getColors <- getColours

#' @param PASS The colour denoting PASS on all plots, in rgb format
#' @param WARN The colour denoting WARN on all plots, in rgb format
#' @param FAIL The colour denoting FAIL on all plots, in rgb format
#' @param MAX The colour denoting the limit of values in rgb format
#' @export
#' @rdname getColours
#' @aliases setColours
setMethod("setColours", "PwfCols", function(object, PASS, WARN, FAIL, MAX){
  new <- object
  if(!missing(PASS)) new@PASS <- PASS
  if(!missing(WARN)) new@WARN <- WARN
  if(!missing(FAIL)) new@FAIL <- FAIL
  if(!missing(MAX)) new@MAX <- MAX

  if(!isValidPwf(new)) {
    warning("Invalid specifications for an object of class PwfCols.\nThe object was not overwritten")
    return(object)
  }

  new
})
#' @export
#' @rdname getColours
#' @aliases setColors
setColors <- setColours

#' @param alpha Numeric(1). Ranges from 0 to 1 by default, but can also be on the range 0 to 255.
#' @export
#' @rdname getColours
#' @aliases setAlpha
setMethod("setAlpha", "PwfCols", function(object, alpha){
  stopifnot(alpha <= 255, alpha >= 0)
  if (alpha > 1) alpha <- alpha/255 # Set to the range [0, 1]
  hexAlpha <- toupper(as.hexmode(floor(alpha*256))) # Convert to hex (upper case)
  hexAlpha <- stringr::str_pad(hexAlpha, width = 2, side = "left", pad = "0")
  setColours(object,
             PASS = paste0(object@PASS, hexAlpha),
             WARN = paste0(object@WARN, hexAlpha),
             FAIL = paste0(object@FAIL, hexAlpha),
             MAX = paste0(object@MAX, hexAlpha))
})

#' @importFrom methods slotNames
#' @export
setMethod(names, "PwfCols", function(x){slotNames(x)})
setMethod("names<-", "PwfCols", function(x){
  warning("The names attribute cannot be set on a PwfCols object")
  x
})

#' The default method for show
#' @param object An object of class PwfCols
setMethod(show, "PwfCols", function(object){
  cat("An object of class PwfCols.\n")
})
