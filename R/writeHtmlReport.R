#' @title Write an HTML Summary Report
#'
#' @description Compiles an HTML report using a supplied template
#'
#' @param fastqcDir A directory containing zipped, or extracted FastQC reports
#' @param template The template file which will be copied into \code{fastqcDir}
#' @param targetsDF A data.frame with at least two columns named \code{Filename} and \code{Label}.
#' The filenames should match the original fastq files, and the labels should be simply alternative
#' labels for these files for convenience.
#' @param overwrite \code{logical}. Overwrite any previous copies of the template file in the destination directory
#' @param quiet \code{logical}. Show or hide markdown output in the Console.
#'
#' @details
#' This will take a user supplied template, or the file supplied with the package and create an HTML
#' summary of all standard FASTQC plots for all files in the supplied directory.
#'
#'
#' @export
writeHtmlReport <- function(fastqcDir, template, targetsDF, overwrite = FALSE, quiet = TRUE){

  # Include checks for the package webshot & PhantomJS
  # Install looks like webshot::install_phantomjs() should work

  if (missing(template)){
    template <- system.file("ngsReports_Fastqc.Rmd", package = "ngsReports")
  }
  if (!file.exists(template)) stop("Could not find template file", template)

  # Copy the template file to fastqcDir
  if (!dir.exists(fastqcDir)) stop("Could not find directory", fastqcDir)
  cp <- file.copy(template, fastqcDir, overwrite = overwrite)
  if(!cp) {
    msg <- paste("Template file was not able to be copied.",
                 "Do you need to set overwrite = TRUE?")
    message(msg)
  }
  file2Knit <- file.path(fastqcDir, basename(template))
  stopifnot(file.exists(file2Knit))

  # Export targets.csv from targets data.frame if supplied
  if (!missing(targetsDF)) {
    chk <- vapply(c("[Ff]ile[Nn]ame", "[Ll]abel"), grepl, logical(1), x = colnames(targetsDF) )
    if (all(chk)) {
      message("Exporting targets.csv")
      readr::write_csv(targetsDF, file.path(fastqcDir, "targets.csv"), append = FALSE)
    }
    else{
      message("Invalid targetsDF. No file will be exported.")
    }
  }

  # Compile the document in the directory
  htmlOut <- file.path(fastqcDir, gsub(".Rmd$", ".html", basename(template)))
  message(paste("Generating", htmlOut, "from template..."))
  rmarkdown::render(file2Knit, output_format = "html_document",
                    output_file = htmlOut, output_dir = dirname(file2Knit),
                    knit_root_dir = dirname(file2Knit), envir = new.env(), quiet = quiet)
  message("done")

  # Return the message, along with invisible(TRUE/FALSE)
  invisible(TRUE)

}
