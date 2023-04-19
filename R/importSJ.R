#' @title Import STAR Splice Junctions
#'
#' @description Import the SJ.out.tab files produced by STAR
#'
#' @details
#' Imports one or more splice-junction output files as produced by STAR.
#' If all are located in separated directories with identical names, be sure
#' to set the argument stripPaths = FALSE
#'
#' All co-ordinates are 1-based, in keeping with the STAR manual
#'
#' @return
#' A `tibble`
#'
#' @param x vector of file paths to SJ.out.tab files
#' @param stripPaths logical(1) Remove directory prefixes from the file paths
#' in x
#'
#' @author Stephen Pederson <stephen.pederson@@adelaide.edu.au>
#'
#' @examples
#' sjFiles <- system.file("extdata", "SJ.out.tab", package = "ngsReports")
#' # Import leaving the complete file path in the column Filename
#' # The argument srtipPaths is set as TRUE by default
#' df <- importSJ(sjFiles, stripPaths = FALSE)
#'
#' @export
importSJ <- function(x, stripPaths = TRUE){

  x <- unique(x) # Remove any  duplicates
  stopifnot(length(x) > 0) # Fail on empty vector
  stopifnot(file.exists(x)) # Check they all exist

  ## Load the data
  data <- suppressWarnings(lapply(x, readLines))
  if (stripPaths) x <- basename(x)
  names(data) <- x
  data <- lapply(data, .splitByTab, firstRowToNames = FALSE)

  ## Check the structure
  allValid <- vapply(data, .isValidSJ, logical(1))
  stopifnot(allValid)
  sjCols <- c(
    "chromosome", "start", "end", "strand", "intron_motif",
    "annotated", "uniquely_mapping", "multi_mapping", "max_overhang"
  )

  ## Name the columns & set as the correct type
  df <- lapply(
    names(data),
    function(x){
      df <- data[[x]]
      colnames(df) <- sjCols
      df[["Filename"]] <- x
      df
    }
  )
  df <- as_tibble(bind_rows(df))
  ## Format all the required columns
  num_cols <- c("start", "end")
  int_cols <- c(
    "strand", "intron_motif", "annotated", "uniquely_mapping", "multi_mapping",
    "max_overhang"
  )
  df[num_cols] <- lapply(df[num_cols], as.numeric)
  df[int_cols] <- lapply(df[int_cols], as.integer)
  df[["annotated"]] <- as.logical(df[["annotated"]])
  ## Tidy up the strands
  strands <- c("*", "+", "-")
  df[["strand"]] <- strands[as.integer(df[["strand"]]) + 1]
  df[["strand"]] <- factor(df[["strand"]], levels = strands)
  ## Recode the motifs
  motifs <-  c(
    "non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT"
  )
  df[["intron_motif"]] <- motifs[df[["intron_motif"]] + 1]
  df[["intron_motif"]] <- factor(df[["intron_motif"]], levels = motifs)
  ## Setup for easy conversion to GRanges
  df[["intron"]] <- paste0(
    df[["chromosome"]], ":",
    df[["start"]], "-", df[["end"]], ":",
    df[["strand"]]
  )

  ## Return the object after creating a dummy variable for R CMD check
  intron <- c()
  dplyr::select(df, Filename, intron, everything())

}

.isValidSJ <- function(df){
  ## Format specified at
  ## https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
  chkColNums <- ncol(df) == 9
  if (!chkColNums) return(FALSE)
  chkStrand <- all(df[[4]] %in% 0:2)
  chkMotif <- all(df[[5]] %in% 0:6)
  chkAnnot <- all(df[[6]] %in% 0:1)
  chkUnique <- sum(is.na(as.numeric(df[[7]]))) == 0
  chkMulti <- sum(is.na(as.numeric(df[[8]]))) == 0
  chkOverhang <- sum(is.na(as.numeric(df[[9]]))) == 0
  all(
    c(
      chkStrand,
      chkMotif,
      chkAnnot,
      chkUnique,
      chkMulti,
      chkOverhang
    )
  )
}
