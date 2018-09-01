#' @title Split elements of a vector into a data.frame
#'
#' @description Split elements of a character vector by the tab separator
#'
#' @details This will split a vector into a data.frame checking that every line
#' has the same number of separators.
#' By default the first element will be set as the column names.
#'
#' This is designed to take input from `readLines()`
#'
#' @param x A character vector
#' @param firstRowToNames logical Should the first element be used for column names
#' @param tab character The string used torepresent the tab symbol
#'
#' @examples
#' x <- c("ColA\tColB", "Value1\tValue2")
#' ngsReports:::splitByTab(x, firstRowToNames = TRUE)
#' ngsReports:::splitByTab(x, firstRowToNames = FALSE)
#'
#' @keywords internal
#'
splitByTab <- function(x, firstRowToNames = TRUE, tab = "\\t"){
    
    stopifnot(is.character(x))
    
    # Check for the tab marker in every line
    linesWithTab <- stringr::str_detect(x, tab)
    if (sum(linesWithTab)!= length(x)) stop("Some elements of x are missing the tab separator")
    
    # Take the first element as defining the number of columns
    nCol <- stringr::str_count(x[1], pattern = tab) + 1
    
    # Count the number of tabs in each line
    nTabs <- stringr::str_count(string = x, pattern = tab)
    if (any(nTabs != (nCol - 1))) stop("Differing number of delimiters in some rows")
    
    if (firstRowToNames){
        
        # Get the first element as a vector of names
        nm <- stringr::str_split_fixed(x[1], pattern = tab, n = nCol)
        
        # Split the remainder
        df <- stringr::str_split_fixed(x[-1], pattern = tab, n = nCol)
        colnames(df) <- nm
    }
    else {
        df <- stringr::str_split_fixed(x, pattern = tab, n = nCol)
    }
    # Return a generic data.frame
    # This leaves prettying to each module
    as.data.frame(df, stringsAsFactors = FALSE)
}
