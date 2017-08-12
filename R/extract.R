#' @title Extract Elements
#'
#' @description Extract elements from FastqcFileList Objects
#'
#' @details Extract elements in a consistent manner with R conventions
#'
#' @param x A FastqcFileList
#' @param i character, logical or integer vector
#' @param j not used
#' @param ... not used
#' @param drop not used
#'
#' @include FastqcFileList.R
#' @include AllGenerics.R
#'
#' @name [
#' @aliases [,FastqcFileList,ANY,missing-method [,FastqcFileList,ANY,missing,ANY-method
#' @rdname extract-methods
#' @export
setMethod("[", signature = c(x = "FastqcFileList", i = "ANY", j = "missing"),
          function(x, i, j, ..., drop){
            if (is.numeric(i) && max(i) > length(x)) stop("Error: subscript out of bounds")
            if (is.logical(i) && length(i) != length(x)) stop("Invalid vector length for subsetting")
            FastqcFileList(x@.Data[i])}
)

#' @name [
#' @aliases [,FastqcFileList,character,missing-method [,FastqcFileList,character,missing,ANY-method
#' @rdname extract-methods
#' @export
setMethod("[", signature = c(x = "FastqcFileList", i = "character", j = "missing"),
          function(x, i, j, ..., drop){
            sub <- match(i, fileName(x))
            if (any(is.na(sub))) {
              msg <- paste("Names not found:", i[is.na(sub)], sep = "\n")
              stop(msg)
            }
            FastqcFileList(x@.Data[sub])}
)

#' @name [[
#' @aliases [[,FastqcFileList,character,missing-method
#' @rdname extract-methods
#' @export
setMethod("[[", signature = c(x = "FastqcFileList", i = "character", j = "missing"),
          function(x, i, j, ...){
            sub <- match(i, fileName(x))
            if (any(is.na(sub))) {
              msg <- paste("Names not found:", i[is.na(sub)], sep = "\n")
              stop(msg)
            }
            x@.Data[[sub]]
          }
)

