#' @title Write an HTML Summary Report
#'
#' @description Compiles an HTML report using a supplied template
#'
#' @param fastqcDir A directory containing zipped, or extracted FastQC reports
#' @param template The template file which will be copied into \code{fastqcDir}
#' @param usePlotly Generate interactive plots?
#' @param species Species/closely related species of sequenced samples
#' @param gcType Is the data "Transcriptomic" or "Genomic" in nature?
#' @param nOver The maximum number of Overrepresented Sequences to show
#' @param nKmer The maximum number of Kmers to show
#' @param targetsDF A data.frame with at least two columns named
#' \code{Filename}  and \code{Label}.
#' The filenames should match the original fastq files, and the labels should
#' be simply alternative labels for these files for convenience.
#' @param overwrite \code{logical}. Overwrite any previous copies of the
#' template file in the destination directory
#' @param quiet \code{logical}. Show or hide markdown output in the Console.
#'
#' @details
#' This will take a user supplied template, or the file supplied with the
#' package and create an HTML summary of all standard FASTQC plots for all
#' files in the supplied directory.
#'
#' @return
#' Silently returns \code{TRUE} and will output a compiled HTML file from the
#' supplied Rmarkdown template file
#'
#' @examples
#' \dontrun{
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names= TRUE)
#' # Copy these files to tempdir() to avoid overwriting
#' # any files in the package directory
#' file.copy(fileList, tempdir())
#' writeHtmlReport(tempdir())
#' }
#'
#' @export
writeHtmlReport <- function(
    fastqcDir, template, usePlotly = TRUE, species = "Hsapiens",
    gcType = c("Genome", "Transcriptome"), nOver = 30, nKmer = 30, targetsDF,
    overwrite = FALSE, quiet = TRUE){

    ## Maybe include checks for the package webshot & PhantomJS
    ## Install looks like webshot::install_phantomjs() should work

    if (missing(template))
        template <- system.file(
            "extdata",
            "ngsReports_Fastqc.Rmd",
            package = "ngsReports"
        )
    if (!file.exists(template))
        stop("Could not find template file", template)
    if (!grepl("Rmd$", template))
        stop("Supplied template must be an rmarkdown file")

    ## Copy the template file to fastqcDir
    if (!dir.exists(fastqcDir)) stop("Could not find directory", fastqcDir)
    cp <- file.copy(template, fastqcDir, overwrite = overwrite)
    if (!cp) {
        msg <- paste(
            "Template file was not able to be copied.",
            "Do you need to set overwrite = TRUE",
            "or check write permissions?"
        )
        stop(msg)
    }
    file2Knit <- file.path(fastqcDir, basename(template))
    stopifnot(file.exists(file2Knit))

    ## Export targets.csv from targets data.frame if supplied
    if (!missing(targetsDF)) {
        chk <- vapply(
            c("[Ff]ile[Nn]ame", "[Ll]abel"),
            FUN = grepl,
            FUN.VALUE = logical(1),
            x = colnames(targetsDF)
        )
        if (all(chk)) {
            message("Exporting targets.csv")
            readr::write_csv(
                targetsDF, file.path(fastqcDir, "targets.csv"), append = FALSE
            )
        }
        else{
            message("Invalid targetsDF. No file will be exported.")
        }
    }

    ## Check the remaining arguments
    stopifnot(is.logical(usePlotly))
    nOver <- suppressWarnings(as.integer(nOver[1]))
    nKmer <- suppressWarnings(as.integer(nKmer[1]))
    stopifnot(!is.na(c(nOver, nKmer)))

    ## Include the namesapce for gcTheoretical as that allows running without
    ## loading the package, which is a pretty common use case
    gcType <- match.arg(gcType)
    avail <- gcAvail(ngsReports::gcTheoretical, gcType)
    species <- match.arg(species, avail$Name)

    ## Compile the document in the directory
    htmlOut <- file.path(fastqcDir, gsub(".Rmd$", ".html", basename(template)))
    message(paste("Generating", htmlOut, "from template..."))
    rmarkdown::render(
        file2Knit,
        output_format = "html_document",
        output_file = basename(htmlOut),
        output_dir = dirname(file2Knit),
        envir = new.env(),
        quiet = quiet,
        params = list(
            usePlotly = usePlotly,
            gcType = gcType,
            species = species,
            nOver = nOver,
            nKmer = nKmer)
    )
    message("done")

    ## Return the message, along with invisible(TRUE/FALSE)
    invisible(TRUE)
}
