#' Get or set the colours in a PwfCols object
#'
#' @description Get or set the colours in a PwfCols object.
#' \code{getColors()} and \code{setColors()} are the same functions as with the non-US spelling
#'
#' @param object An object of class PwfCols
#' @param PASS The colour used to signify a PASS
#' @param WARN The colour used to signify a WARN
#' @param FAIL The colour used to signify a FAIL
#' @param MAX The colour used as an extreme FAIL
#' @param alpha Either a value between 0 and 1, or in the range [0, 255].
#' See \code{\link{col2rgb}}
#'
#' @details All colours must be specified in rgb format.
#'
#' Calling the function \code{plot()} on an object of class \code{PwfCols} will draw a pie chart
#' showing the assigned colours.
#'
#' @return An object of class \code{PwfCols}
#'
#' @export
#' @rdname getColours
setMethod("getColours", "PwfCols", function(object){
  vals <- c(object@PASS,
            object@WARN,
            object@FAIL,
            object@MAX)
  names(vals) <- c("PASS", "WARN", "FAIL", "MAX")
  vals
})

#' @export
#' @rdname getColours
setMethod("getColors", "PwfCols", function(object){
  vals <- c(object@PASS,
            object@WARN,
            object@FAIL,
            object@MAX)
  names(vals) <- c("PASS", "WARN", "FAIL", "MAX")
  vals
})

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
#' @aliases setColours
setMethod("setColors", "PwfCols", function(object, PASS, WARN, FAIL, MAX){
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

#' @importFrom stringr str_pad
#' @export
#' @rdname getColours
#' @aliases setAlpha,PwfCols-method
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
