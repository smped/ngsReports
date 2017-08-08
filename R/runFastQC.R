#' A wrapper for the bash shell command fastqc.
#'
#' @description Takes a FastqFileList
#'
#' @details This is a simple wrapper function for controlling & running \code{fastqc} from within R.
#' This can be very useful for controlling & documenting an entire pipeline from within knitr to produce a simple report
#'
#' Only the common functinoality of FastQC is implemented,
#' for more fine detail pease call FastQC directly.
#'
#' @param object A \code{\link{FastqFileList}}
#' @param outPath The path to write the FastQC reports. Must exist as it will not be created when calling this function
#' @param exec The location of the fastqc executable.
#' @param threads The number of threads to run in parallel
#' @param casava logical. Sets the \code{--casava} flag
#' @param nofilter logical. Sets the \code{--nofilter} flag
#' @param extract logical. Extract the zip files on completion of the report
#' @param nogroup logical. Sets the grouping of bases for reads longer than 50bp
#' @param contaminants Path to an alternate file with contaminants.
#' The structure of the file will not be checked.
#' Refer to the \code{fastqc} help page for more details
#' @param adapters Path toa file listing adapters to search for.
#' The structure of the file will not be checked.
#' Refer to the \code{fastqc} help page for more details
#' @param kmers An integer between 2 and 10
#' @param overWrite logical. Overwrite any existing reports
#'
#' @return An object of class \code{FastqcFileList}
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @seealso \code{\link{FastqFileList}}
#'
#' @importFrom parallel detectCores
#'
#' @export
runFastQC <- function(object, outPath, exec = "/usr/local/bin/fastqc",
                      threads=1L, casava = FALSE, nofilter = FALSE, extract = FALSE,
                      nogroup = FALSE, contaminants = c(), adapters = c(), kmers = 7,
                      overWrite = FALSE){

  # Check the FastqFileList
  stopifnot(class(object) == "FastqFileList")
  nm <- names(object)
  if (any(!grepl("(fastq|fastq.gz|fq|fq.gz)$", Rsamtools::path(object)))) stop("Files can only contain the fasqt|fq suffix")
  # Make sure they are all of the same format
  suffix <- gsub(".+(fastq|fastq.gz|fq|fq.gz)$", "\\1", nm)
  suffix <- paste0(".", unique(suffix))
  stopifnot(length(suffix) == 1)

  # Check the output path exists
  stopifnot(file.exists(outPath))
  if (length(outPath) > 1){
    outPath <- outPath[1]
    message("Multiple output paths specified.\nOnly the first will be used.")
  }

  # Check the executable exists
  stopifnot(file.exists(exec))

  # Set the expected output names
  fqcNames <- file.path(outPath, gsub(suffix, "_fastqc", nm))
  if (!extract) fqcNames <- paste0(fqcNames, ".zip")
  if (!overWrite && all(file.exists(fqcNames))) {
    message("All reports exist and were not overwritten")
    return(FastqcFileList(fqcNames))
  }

  # It might be worth addiing a section here to not overwrite any existing files,
  # and to run fastQC on the remaining files
  # The entire set would need to be kept for the final output though

  # Set the arguments to the function call
  maxCores <- parallel::detectCores()
  threads <- paste("-t", min(threads, maxCores))
  stopifnot(is.logical(casava), is.logical(nogroup))
  if (casava){
    casava <- dplyr::if_else(nofilter, "--casava --nofilter", "--casava")
  }
  else {
    casava <- c()
  }
  extract <- dplyr::if_else(extract, "--extract", "--noextract")
  if (nogroup) {
    warning("Setting the option '--nogroup' may cause FastQC to become unstable")
    nogroup <- "--nogroup"
  }
  else{
    nogroup <- c()
  }
  if (!is.null(contaminants)) {
    stopifnot(file.exists(contaminants))
    contaminants <- paste("-c", contaminants)
  }
  if (!is.null(adapters)){
    stopifnot(file.exists(adapters))
    adapters <- paste("-a", adapters)
  }
  kmers <- as.integer(kmers)[1]
  stopifnot(!is.na(kmers), kmers %in% 2:10)
  kmers <- paste("-k", kmers)

  args <- paste("-o", outPath, threads, casava, extract, nogroup, contaminants, adapters, kmers)
  args <- gsub(" +", " ", args) #Remove any double spaces
  files <- paste(Rsamtools::path(object), collapse = " ")

  # Run the command
  message("Executing the command '", paste(exec, args), "'")
  system2(exec, paste(args, files))

  # Check that everything was written as expected
  if (any(!file.exists(fqcNames))) warning("Reports not written for:\n", fqcNames[!file.exists(fqcNames)])

  # Get the version number for the output
  v <- system2(exec, "-v", stdout = TRUE)

  # Return a FastqcFileList, until extracting the summary has been incorporated into the
  # getFastqcData methods for a FastqcFileList
  FastqcFileList(fqcNames)

}
