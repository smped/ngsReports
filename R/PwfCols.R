#' The PwfCols class and associated methods
#'
#' @description Define the PwfCols class and associated methods
#'
#' @details
#' This is an object of with four colours in components named PASS, WARN, FAIL
#' and MAX.
#' Used to indicate these categories as defined on the standard plots from
#' fastqc.
#'
#' @slot PASS A vector of length 1, defining the colour for PASS in rgb format.
#' Defaults to rgb(0, 0.8, 0)
#' @slot WARN A vector of length 1, defining the colour for WARN in rgb format.
#' Defaults to rgb(0.9, 0.9, 0.2)
#' @slot FAIL A vector of length 1, defining the colour for FAIL in rgb format.
#' Defaults to rgb(0.8, 0.2, 0.2)
#' @slot MAX A vector of length 1, defining the colour for an extreme FAIL or
#' NA in rgb format. Defaults to rgb(1, 1, 1)
#'
#' @include validationFunctions.R
#' @aliases PwfCols
setClass("PwfCols", slots = c(
    PASS = "character",
    WARN = "character",
    FAIL = "character",
    MAX = "character")
)
setValidity("PwfCols", .isValidPwf)

# The default method for show. Doesn't need exporting
setMethod("show", "PwfCols",  function(object){
    cat("An object of class PwfCols.\n")
})