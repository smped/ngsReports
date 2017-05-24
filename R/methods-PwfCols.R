#' Define methods for objects of class PwfCols
#'
#'
#' @include AllClasses.R
#' @include AllGenerics.R
#'
#' @docType methods
#'
#' @importFrom graphics plot
#'
#' @export
#' @rdname PwfCols-methods
setMethod(names, "PwfCols", function(x){slotNames(x)})

setMethod(show, "PwfCols", function(object){
  cat("An object of class PwfCols.\n")
  })

#' @export
#' @rdname PwfCols-methods
setMethod("plot", "PwfCols", function(x, ...){
  pie(rep(1,4), labels = names(x), col = getColours(x), ...)
})
