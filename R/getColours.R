#' Get or set the colours in a PwfCols object
#'
#' @description Get or set the colours in a PwfCols object
#'
#' @param object An object of class PwfCols
#' @param PASS The colour used to signify a PASS
#' @param WARN The colour used to signify a WARN
#' @param FAIL The colour used to signify a FAIL
#' @param MAX The colour used as an extreme FAIL
#'
#' @details All colours must be specified in rgb format.
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
