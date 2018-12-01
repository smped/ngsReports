#' @title Work with objects of class PwfCols
#'
#' @description Get and modify colours from objects of class PwfCols
#'
#' @details
#' Use \code{getColours} to obtain the colours in an object of class PwfCols.
#'
#' These can be modified using the functions \code{setColours} and
#' \code{setAlpha}
#'
#' @param object An object of class PwfCols
#'
#' @return
#' getColours will return a character vector of colours coresponding to
#'  PASS/WARN/FAIL
#'
#' @examples
#' getColours(pwf)
#'
#' # How to add transparency
#' pwf2 <- setAlpha(pwf, 0.1)
#' getColours(pwf2)
#'
#' @include AllGenerics.R
#' @importFrom methods slotNames
#' @name getColours
#' @aliases getColours,PwfCols-method
#' @rdname getColours-methods
#' @export
setMethod("getColours", "PwfCols", function(object){
    vals <- c(
        object@PASS,
        object@WARN,
        object@FAIL,
        object@MAX
    )
    names(vals) <- c("PASS", "WARN", "FAIL", "MAX")
    vals
})

#' @param PASS The colour denoting PASS on all plots, in rgb format
#' @param WARN The colour denoting WARN on all plots, in rgb format
#' @param FAIL The colour denoting FAIL on all plots, in rgb format
#' @param MAX The colour denoting the limit of values in rgb format
#
#' @return
#' setColours will return an object of class PwfCols
#'
#' @importFrom methods slotNames
#' @name setColours
#' @aliases setColours,PwfCols-method setColours
#' @rdname getColours-methods
#' @export
setMethod("setColours", "PwfCols", function(object, PASS, WARN, FAIL, MAX){
    new <- object
    if (!missing(PASS)) new@PASS <- PASS
    if (!missing(WARN)) new@WARN <- WARN
    if (!missing(FAIL)) new@FAIL <- FAIL
    if (!missing(MAX)) new@MAX <- MAX

    if (!.isValidPwf(new)) {
        wn <- paste(
            "Invalid specifications for an object of class PwfCols.",
            "The object was not overwritten", sep = "\n"
        )
        warning(wn)
        return(object)
    }

    new
})

#' @name setAlpha
#' @param alpha Numeric(1). Ranges from 0 to 1 by default, but can also be on
#' the range 0 to 255.
#' @return
#' setAlpha will return an object of class PwfCols
#' @importFrom methods slotNames
#' @export
#' @rdname getColours-methods
#' @aliases setAlpha,PwfCols-method
setMethod("setAlpha", "PwfCols", function(object, alpha){
    stopifnot(alpha <= 255, alpha >= 0)

    ## Get the alpha value in hex
    if (alpha > 1) alpha <- alpha/255 # Set to the range [0, 1]
    hexAlpha <- toupper(as.hexmode(floor(alpha*256))) # Convert to hex
    hexAlpha <- stringr::str_pad(hexAlpha, width = 2, side = "left", pad = "0")

    ## Get the colours without any existing alpha
    oldCols <- getColours(object)
    oldCols <- substr(oldCols, start = 1, stop = 7)
    newCols <- paste0(oldCols, hexAlpha)
    names(newCols) <- names(oldCols)

    ## Now assign them to a new object
    args <- c(Class = "PwfCols", as.list(newCols))
    do.call(new, args)
})
