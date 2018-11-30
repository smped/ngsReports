#' @title Calculate GC content from a fasta file
#' 
#' @description Calculate GC content from a set of DNA sequences provided by a 
#' fasta file
#' 
#' @details Used to obtain theoretical GC content for a given 
#' genome/transcriptome
#' 
#' @param Fastafile A fasta file containing a genome or transcriptome of 
#' interest
#' @param n the number of fragments to simulate
#' @param bp the fragment length used in estimation
#' 
#' @return A data.frame with the GC percentage and the exected frequency
#' 
#' @author Hien To <totuhien@gmail.com>
#' 
#' @keywords internal
gcFromFasta <- function(Fastafile, n=1e+6, bp=100){
    # Load in the reference sequences
    stopifnot(file.exists(Fastafile))
    ref <- readDNAStringSet(filepath = Fastafile)
    # Define a temporary export file for generateDistn
    file.gc <- file.path(tempdir(), "gc.txt")
    fastqcTheoreticalGC::generateDistn(ref, file = file.gc, n = n, bp = bp)
    # Load in the results & delete the tempfile
    gc <- read.table(file.gc, header = FALSE, stringsAsFactors = FALSE)
    gc <- data.frame(gc)
    file.remove(file.gc)
    # Format the output
    colnames(gc) <- c("GC_Content","Freq")
    gc$Freq <- gc$Freq/sum(gc$Freq)
    gc
}