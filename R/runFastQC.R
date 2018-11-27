#' A wrapper for the bash shell command fastqc.
#'
#' @description Takes a FastqFileList
#'
#' @details This is a simple wrapper function for controlling & running 
#' \code{fastqc} from within R.
#' This can be very useful for controlling & documenting an entire pipeline from 
#' within knitr to produce a simple report
#'
#' Only the common functionality of FastQC is implemented,
#' for more fine detail pease call FastQC directly.
#'
#' @param object A \code{\link{FastqFileList}}
#' @param outPath The path to write the FastQC reports. Must exist as (for safety)
#' it will not be created when calling this function
#' @param exec The location of the fastqc executable.
#' @param threads The number of threads to run in parallel
#' @param casava logical. Sets the \code{--casava} flag
#' @param nofilter logical. Sets the \code{--nofilter} flag
#' @param extract logical. Extract the zip files on completion of the report
#' @param nogroup logical. Sets the grouping of bases for reads longer than 50bp
#' @param min_length integer. Sets an artificial lower limit on the length of the 
#' sequence to be shown in the report.  
#' @param contaminants Path to an alternate file with contaminants.
#' The structure of the file will not be checked.
#' Refer to the \code{fastqc} help page for more details
#' @param adapters Path toa file listing adapters to search for.
#' The structure of the file will not be checked.
#' Refer to the \code{fastqc} help page for more details
#' @param kmers An integer between 2 and 10
#'
#' @return An object of class \code{FastqcFileList}
#'
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @seealso \code{\link{FastqFileList}}
#' 
#' @examples 
#' \dontrun{
#' library(ShortRead)
#' sp <- SolexaPath(system.file('extdata', package='ShortRead'))
#' fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
#' f <- FastqFile(fl)
#' # This requires a working installation of FastQC
#' fqcFile <- runFastQC(f, outPath = tempdir())
#' }
#'
#' @importFrom parallel detectCores
#' @importClassesFrom ShortRead FastqFile FastqFileList
#' @importClassesFrom Rsamtools BamFile BamFileList
#'
#' @export
#' @rdname runFastQC-methods
setGeneric("runFastQC", function(object, outPath, threads=1L, casava = FALSE, 
                                 nofilter = FALSE, extract = FALSE, 
                                 nogroup = FALSE, min_length = 1, 
                                 contaminants = c(), adapters = c(), kmers = 7,
                                 exec){
    standardGeneric("runFastQC")
})
#' @aliases runFastQC,character-method
#' @rdname runFastQC-methods
#' @export
setMethod("runFastQC", "character",
          function(object, outPath, threads=1L, casava = FALSE, 
                   nofilter = FALSE, extract = FALSE, 
                   nogroup = FALSE, min_length = 1, 
                   contaminants = c(), adapters = c(), kmers = 7,
                   exec){
              
             fq <- all(grepl("(fq|fastq|fq.gz|fastq.gz|txt)$", basename(object)))
             bam <- all(grepl("bam$", basename(object)))
             if (fq) {
                 if (length(object) == 1) {
                     object <- ShortRead::FastqFile(object)
                 }
                 else {
                     object <- ShortRead::FastqFileList(object)
                 }
             }
             if (bam) {
                 if (length(object) == 1) {
                     object <- Rsamtools::BamFile(object)
                 }
                 else{
                     object <- Rsamtools::BamFileList(object)
                 }
             }
             if (!bam & !fq) {
                 warning("File format could not be identified. FastQC will not be run")
                 return(NULL)
             }
             runFastQC(object, outPath, threads, casava, nofilter, 
                       extract, nogroup, min_length, contaminants = c(), 
                       adapters, kmers, exec)
          })
#' @aliases runFastQC,FastqFile-method
#' @rdname runFastQC-methods
#' @export
setMethod("runFastQC", "FastqFile",
          function(object, outPath, threads=1L, casava = FALSE, 
                   nofilter = FALSE, extract = FALSE, 
                   nogroup = FALSE, min_length = 1, 
                   contaminants = c(), adapters = c(), kmers = 7,
                   exec){
              
              runFastQC.default(object, outPath, threads, casava, nofilter, 
                                extract, nogroup, min_length, contaminants = c(), 
                                adapters, kmers, exec, fileType = "FastqFile") 
          }
)
#' @aliases runFastQC,FastqFileList-method
#' @rdname runFastQC-methods
#' @export
setMethod("runFastQC", "FastqFileList",
          function(object, outPath, threads=1L, casava = FALSE, 
                   nofilter = FALSE, extract = FALSE, 
                   nogroup = FALSE, min_length = 1, 
                   contaminants = c(), adapters = c(), kmers = 7,
                   exec){
        
              runFastQC.default(object, outPath, threads, casava, nofilter, 
                                extract, nogroup, min_length, contaminants = c(), 
                                adapters, kmers, exec, fileType = "FastqFileList")      
              
          }
)
#' @aliases runFastQC,BamFile-method
#' @rdname runFastQC-methods
#' @export
setMethod("runFastQC", "BamFile",
          function(object, outPath, threads=1L, casava = FALSE, 
                   nofilter = FALSE, extract = FALSE, 
                   nogroup = FALSE, min_length = 1, 
                   contaminants = c(), adapters = c(), kmers = 7,
                   exec){
              
              runFastQC.default(object, outPath, threads, casava, nofilter, 
                                extract, nogroup, min_length, contaminants = c(), 
                                adapters, kmers, exec, fileType = "BamFile") 
          }
)
#' @aliases runFastQC,BamFileList-method
#' @rdname runFastQC-methods
#' @export
setMethod("runFastQC", "BamFileList",
          function(object, outPath, threads=1L, casava = FALSE, 
                   nofilter = FALSE, extract = FALSE, 
                   nogroup = FALSE, min_length = 1, 
                   contaminants = c(), adapters = c(), kmers = 7,
                   exec){
              
              runFastQC.default(object, outPath, threads, casava, nofilter, 
                                extract, nogroup, min_length, contaminants = c(), 
                                adapters, kmers, exec, fileType = "BamFileList") 
          }
)

runFastQC.default <- function(object, outPath, threads, casava, nofilter, 
                              extract, nogroup, min_length, contaminants, 
                              adapters, kmers, exec, fileType){
    
    stopifnot(!missing(fileType))
    
    # Check the output path exists
    stopifnot(file.exists(outPath))
    if (length(outPath) > 1){
        outPath <- outPath[1]
        message("Multiple output paths specified.\nOnly the first will be used.")
    }
    
    # Check the executable exists.
    if (missing(exec)) exec <- Sys.which("fastqc")
    if (file.exists(exec)){
        message("Found FastQC executable: ", exec)
    }
    else {
        stop("FastQC executable not found")
    }
    
    # Get the version
    v <- system2(exec, "-v", stdout = TRUE)
    v <- gsub(".+0.11.(.+)", "\\1", v)
    if (as.integer(v) > 5){
        min_length <- paste("--min_length", min_length)
    }
    else{
        message("min_length cannot be set for FastQC < 0.11.6")
        min_length <- c()
    }
    
    # Set the arguments to the function call
    maxCores <- parallel::detectCores() - 1
    if (maxCores < threads) message("Too many threads. Resetting to", maxCores)
    threads <- paste("-t", min(threads, maxCores))
    stopifnot(is.logical(casava), is.logical(nogroup))
    if (casava){
        casava <- ifelse(nofilter, "--casava --nofilter", "--casava")
    }
    else {
        casava <- c()
    }
    extract <- dplyr::if_else(extract, "--extract", "--noextract")
    if (nogroup) {
        message("Setting the option '--nogroup' may cause FastQC to become unstable")
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
    
    args <- paste("-o", outPath, threads, casava, extract, nogroup, 
                  contaminants, adapters, kmers)
    args <- gsub(" +", " ", args) #Remove any double spaces
    files <- paste(BiocGenerics::path(object), collapse = " ")
    
    # Run the command
    message("Executing the command '", paste(exec, args), "'")
    system2(exec, paste(args, files))
    
    # Get the files and make sure there are no html files returned
    fqcNames <- list.files(outPath, pattern = "fastqc", 
                           full.names = TRUE)
    fqcNames <- grep(fqcNames, 
                     pattern = "html", invert = TRUE, value = TRUE)
    
    
    # Now define the format as required
    if (fileType %in% c("FastqFileList", "BamFileList")){
        # Return a FastqcFileList
        ffl <- FastqcFileList(fqcNames)
        # Now return the files that have been run in the same order
        m <- pmatch(gsub("(.+)\\..+", "\\1", names(object)), fileName(ffl))
        if (anyNA(m)) {
            warning("NA values indicate Fastqc files may be missing for one or more samples")
            m <- m[!is.na(m)]
        }
        out <- ffl[m]    
    }
    if (fileType %in% c("FastqFile", "BamFile")){
        # Return a FastqcFileList
        ffl <- FastqcFileList(fqcNames)
        # Now return the file that has been run 
        nm <- basename(path(object))
        m <- pmatch(gsub("(.+)\\..+", "\\1", nm), fileName(ffl))
        out <- ffl[[m]]    
    }
    
    # Close any open connections. Nope. This spits errors
    #if (isOpen(object)) close(object)
    out
}
