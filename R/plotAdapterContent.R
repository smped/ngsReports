#' @title Draw an Adapter Content Plot
#'
#' @description Draw an Adapter Content Plot across one or more FASTQC reports
#'
#' @details
#' This extracts the Adapter_Content from the supplied object and generates a ggplot2 object,
#' with a set of minimal defaults.
#' The output of this function can be further modified using the standard ggplot2 methods.
#'
#' @param x Can be a \code{FastqcFile}, \code{FastqcFileList}, \code{FastqcData},
#' \code{FastqcDataList} or path
#' @param facet logical. Determines whether to use \code{facet_wrap()} or not
#' @param nc Number of columns to use when faceting by Adapter Type.
#' Must be supplied if using facets.
#' @param ylim A vector providing limits for the y-axis. Defaults to \code{ylim = c(0, 100)}
#' @param th A ggplot2 theme object. Defaults to \code{theme_bw()}
#' @param trimNames logical. Remove the file suffix from the legend.
#' @param ... Not used
#'
#' @return A standard ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#'
#' @export
plotAdapterContent <- function(x, facet = TRUE, nc, ylim, th, trimNames = TRUE, ...){

  if (missing(ylim)) {
    message("The 'ylim' argument was not specified. The default of (0, 100) will be applied.")
    ylim <- c(0, 100)
  }
  else {
    if (any(length(ylim) != 2, !is.numeric(ylim))) stop("The argument ylim must be a numeric vector of length 2")
  }
  stopifnot(is.logical(facet))

  # Get the AdapterContent
  ac <- Adapter_Content(x)

  # # Check the alternate names
  # if (!missing(altNames)) {
  #   if (any(!ac$Filename %in% names(altNames))) {
  #     warning("Supplied vector of alternate names does not match and will be ignored.")
  #   }
  #   else{
  #     ac$Filename <- altNames[ac$Filename]
  #   }
  # }
  if (trimNames) ac$Filename <- gsub("(.+)\\.(fastq|fq).*", "\\1", ac$Filename)

  # Find the average of the positions
  # This saves problems with non-numeric ordering
  # as well as dealing with potential conflicts in fastqc files themselves
  ac <- dplyr::mutate(ac,
                      Start = gsub("([0-9]*)-[0-9]*", "\\1", Position),
                      Start = as.integer(Start),
                      End = gsub("[0-9]*-([0-9]*)", "\\1", Position),
                      End = as.integer(End),
                      Position = 0.5*(Start + End))
  ac <- dplyr::select(ac, -Start, -End)

  # Change to long form and remove the _ symbols between words
  ac <- reshape2::melt(ac, id.vars = c("Filename", "Position"),
                       value.name = "Percent", variable.name = "Type")
  ac <- mutate(ac, Type = gsub("_", " ", Type))

  # Create the basic plot
  acPlot <- ggplot2::ggplot(ac, ggplot2::aes(x = Position, y = Percent, colour = Filename)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::labs(x = "Position in read (bp)")

  # Add the basic customisations
  if (facet) {
    if (missing(nc)){
      message("The argument nc was not supplied. Using the default of nc = 2")
      nc <- 2
    }
    acPlot <- acPlot + ggplot2::facet_wrap(~Type, ncol = nc)
  }
  # Apply theme_bw() if missing, or an invalid theme is supplied
  if (missing(th)) {
    th <- ggplot2::theme_bw()
  }
  else{
    if (any(!c("theme", "gg") %in% class(th))) {
      warning("Theme supplied as th is not a ggplot theme and will be ignored")
      th <- ggplot2::theme_bw()
    }
  }
  acPlot <- acPlot + th

  # And draw the plot
  acPlot
}
