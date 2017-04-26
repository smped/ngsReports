#' @title Colours for PASS/WARN/FAIL
#'
#' @description Default colours for PASS/WARN/FAIL
#'
#' @details \code{pwfCols} must be a vector with four colours in rgb format.
#' Colours correspond to PASS, WARN and FAIL from the FASTQC reports,
#' with the additional colour denoted as X to indicate an eXtreme FAIL.
#' In order, these colours in the defulat vector are green (\code{rgb(0, 0.8,0)}),
#' yellow (\code{rgb(0.9, 0.9, 0.2)}), red (\code{rgb(0.8, 0.2, 0.2)}) and white
#' (\code{rgb(1, 1, 1)})
#'
#' @param x a vector to be tested as a valid pwfCols vector
#' @param pwfCols a default vector supplied with the package. Described in Details below
#'
#' @return validPwf returns a \code{logical} vector of length 1
#'
#'
#' @rdname pwfCols
#' @export
validPwf <- function(x){
  if (length(x) != 4) return(FALSE)
  if (any(substr(x, 1, 1) != rep("#", 4))) return(FALSE)
  if (any(nchar(x) != rep(7, 4))) return(FALSE)
  if (any(grepl("[G-Z]", x))) return(FALSE)

  TRUE
}

#' @rdname pwfCols
#' @export
pwfCols <- c(PASS = rgb(0, 0.8,0), WARN = rgb(0.9, 0.9, 0.2), FAIL = rgb(0.8, 0.2, 0.2),
             X = rgb(1, 1, 1))
