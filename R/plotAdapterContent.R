#' @title Draw an Adapter Content Plot
#'
#' @description Draw an Adapter Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the Adapter_Content from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' To set any faceting parameters outside the function, such as \code{scales = "free_y"},
#' set \code{facet = FALSE} and apply \code{facet_wrap()} after the call to \code{plotAdapterContent}.
#' Note that the only column available for faceting with be \code{Type},
#' so the formula \code{~Type} will still need to be applied.
#'
#' Preset axis limits can also be overwritten easily by adding a call to \code{scale_y_continuous()}
#' after the call to \code{plotAdapterContent}.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param facet logical. Determines whether to use \code{facet_wrap()} or not
#' @param nc Number of columns to use when faceting by Adapter Type.
#' @param ylim A \code{numeric vector} providing limits for the y-axis.
#' Defaults to \code{ylim = c(0, 100)}
#' @param trimNames \code{logical}. Remove the file suffix from the names displyed in the legend.
#' @param pattern \code{character}.
#' Contains a regular expression which will be captured from fileNames.
#' The default will capture all text preceding .fastq/fastq.gz/fq/fq.gz
#' @param type A regular expression used to filter which adapter(s) are plotted
#' @param ... Used to pass additional arguments, such as \code{invert} or \code{ignore.case},
#' to \code{grep} when selecting Adapter \code{type}
#'
#' @return A standard ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom stringr str_detect
#'
#' @export
plotAdapterContent <- function(x, subset,
                               facet = TRUE, nc = 2, ylim = c(0, 100),
                               trimNames = TRUE, pattern = "(.+)\\.(fastq|fq).*", type, ...){

  # A basic cautionary check
  stopifnot(grepl("(Fastqc|character)", class(x)))

  if (missing(subset)){
    subset <- rep(TRUE, length(x))
  }
  stopifnot(is.logical(subset))
  stopifnot(length(subset) == length(x))
  stopifnot(is.logical(trimNames))
  stopifnot(is.numeric(ylim))
  stopifnot(length(ylim) == 2)
  stopifnot(is.logical(facet))

  # Get the AdapterContent
  x <- x[subset]
  df <- tryCatch(Adapter_Content(x))

  # Check the pattern contains a capture
  if (trimNames && stringr::str_detect(pattern, "\\(.+\\)")) {
    df$Filename <- gsub(pattern[1], "\\1", df$Filename)
    # These need to be checked to ensure non-duplicated names
    if (length(unique(df$Filename)) != length(x)) stop("The supplied pattern will result in duplicated filenames, which will not display correctly.")
  }

  # Find the average of the positions
  # This saves problems with non-numeric ordering
  # as well as dealing with potential conflicts in fastqc files themselves
  df <- dplyr::mutate(df,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Position),
                      Start = as.integer(Start),
                      End = gsub("[0-9]*-([0-9]*)", "\\1", Position),
                      End = as.integer(End),
                      Position = 0.5*(Start + End))
  df <- dplyr::select(df, -Start, -End)

  # Change to long form and remove the _ symbols between words
  df <- reshape2::melt(df, id.vars = c("Filename", "Position"),
                       value.name = "Percent", variable.name = "Type")
  df <- mutate(df, Type = gsub("_", " ", Type))

  # Restrict to a given type if requested
  if (!missing(type)) {
    keep <- grep(pattern = type, x = df$Type, ...)
    df <- df[keep,]
    if(nrow(df) == 0) stop("No adapters matching the supplied type were found")
  }

  # Create the basic plot
  acPlot <- ggplot2::ggplot(df, ggplot2::aes(x = Position, y = Percent, colour = Filename)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::labs(x = "Position in read (bp)",
                  y = "Percent (%)") +
    ggplot2::theme_bw()

  # Add the basic customisations
  if (facet) {
    if (!is.numeric(nc)){
      message("The argument nc was not supplied correctly. Using the default of nc = 2")
      nc <- 2
    }
    acPlot <- acPlot + ggplot2::facet_wrap(~Type, ncol = nc[1])
  }

  # And draw the plot
  acPlot

}
