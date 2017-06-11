#' @title Get the maximum Adapter Content
#'
#' @description Get the maximum Adapter Content across one or more FASTQC reports
#'
#' @details This will extract the \code{Adapter_Content} from the supplied object,
#' and provide a \code{data_frame} with the final value for each file
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param digits \code{numeric}.
#' The output will have the percentages rounded to the specified number of digits
#' @param asPercent \code{logical}.
#' Format the values as percentages with the added \code{\%} symbol
#'
#' @return A \code{data_frame} object containing the percent of reads with each adapter
#' type at the final position
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom scales percent
#'
#' @export
maxAdapterContent <- function(x, digits = 2, asPercent = TRUE){

  # Get the AdapterContent
  ac <- tryCatch(Adapter_Content(x))

  # Perform the summary
  ac %<>%
    reshape2::melt(id.vars = c("Filename", "Position"),
                   variable.name = "Type")%>%
    dplyr::group_by(Filename, Type) %>%
    dplyr::summarise(value = max(value))

  # Format the output
  if (!is.numeric(digits)) {
    message("The digits argument was not supplied as a numeric value and will be ignored.")
  }
  else {
    digits <- floor(digits)[1] # Silently ignore any additional values
    ac %<>% dplyr::mutate(value = round(value, digits))
  }
  if (asPercent) ac %<>% dplyr::mutate(value = scales::percent(0.01*value))

  ac %>%
    reshape2::dcast(Filename~Type, value.var = "value") %>%
    tibble::as_tibble()

}
