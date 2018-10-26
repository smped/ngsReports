#' @title Get the maximum Adapter Content
#'
#' @description Get the maximum Adapter Content across one or more FASTQC reports
#'
#' @details This will extract the \code{Adapter_Content} from the supplied object,
#' and provide a \code{tibble} with the final value for each file
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param digits \code{numeric}.
#' The output will have the percentages rounded to the specified number of digits
#' @param asPercent \code{logical}.
#' Format the values as percentages with the added \code{\%} symbol
#'
#' @return A \code{tibble} object containing the percent of reads with each adapter
#' type at the final position
#'
#'
#' @export
maxAdapterContent <- function(x, digits = 2, asPercent = TRUE){
    
    stopifnot(is.numeric(digits) && is.logical(asPercent))
    
    # Get the AdapterContent
    ac <- tryCatch(Adapter_Content(x))
    
    # Perform the summary
    ac <- reshape2::melt(ac, id.vars = c("Filename", "Position"),
                         variable.name = "Type")
    ac <- dplyr::group_by(ac, Filename, Type)
    ac <- dplyr::summarise_at(ac, dplyr::vars("value"), dplyr::funs("max"))
    
    # Format the output
    digits <- floor(digits)[1] # Silently ignore any additional values
    ac$value <- round(ac$value, digits)
    
    if (asPercent) ac$value <- scales::percent(0.01*ac$value)
    
    ac <- reshape2::dcast(ac, Filename~Type, value.var = "value")
    tibble::as_tibble(ac)
    
}
