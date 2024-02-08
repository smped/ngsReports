#' @title Summarise Overrepresented Sequences
#'
#' @description
#' Summarise the Overrepresented sequences found in one or more QC files
#'
#' @details
#' This function prepares a useful summary of all over-represented sequences
#' as reported by either fastp or FastQC
#'
#'
#' @param x An object of a suitable class
#' @param min_count Filter sequences with counts less than this value, both
#' before and after filtering
#' @param step Can be 'Before', 'After' or both to obtain data from the
#' Before_filtering or After_filtering modules
#' @param vals Values to use for creating summaries across multiple files. For
#' FastpDataList objects these can be "count" and/or "rate", whilst for
#' FastqcDataList objects these values can be "Count" and/or "Percentage"
#' @param fn Functions to use when summarising values across multiple files
#' @param pattern Regular expression to filter the Possible_Source column by
#' @param by character vector of columns to summarise by. See [dplyr::summarise]
#' @param ... Not used
#'
#' @return A tibble
#'
#' Tibble columns will vary between Fastp*, FastqcDataList and FastqcData
#' objects. Calling this function on list-type objects will attempt to
#' summarise the presence each over-represented sequence across all files.
#'
#' In particular, FastqcData objects will provide the requested summary
#' statistics across all sequences within a file
#'
#'
#' @examples
#' ## For operations on a FastpData object
#' f <- system.file("extdata/fastp.json.gz", package = "ngsReports")
#' fp <- FastpData(f)
#' summariseOverrep(fp, min_count = 100)
#'
#' ## Applying the function to a FastqcDataList
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#' fdl <- FastqcDataList(fl)
#' summariseOverrep(fdl)
#'
#' # An alternative viewpoint can be obtained using
#' fdl |> lapply(summariseOverrep) |> dplyr::bind_rows()
#'
#'
#' @name summariseOverrep
#' @rdname summariseOverrep-methods
#' @export
setGeneric(
  "summariseOverrep",
  function(x, ...) standardGeneric("summariseOverrep")
)
#' @rdname summariseOverrep-methods
#' @importFrom dplyr bind_rows if_all
#' @importFrom tidyr unnest pivot_wider
#' @importFrom tidyselect all_of starts_with
#' @export
setMethod(
  "summariseOverrep", signature = "FastpData",
  function(x, step = c("Before", "After"), min_count = 0, ...) {
    step <- match.arg(step, several.ok = TRUE)
    modNames <- paste0(step, "_filtering")
    modData <- lapply(modNames, \(i) getModule(x, i))
    names(modData) <- step
    modData <- lapply(modData, bind_rows, .id = "reads")
    df <- bind_rows(modData, .id = "step")
    if (!"overrepresented_sequences" %in% names(df)) stop(
      "No overrepresented sequence data available"
    )
    df <- dplyr::select(
      df, all_of(c("fqName", "step", "reads")), starts_with("overrep")
    )
    df <- unnest(df, starts_with("overrep"))
    df <- pivot_wider(
      df,
      names_from = "step", values_from = c("count", "rate"), values_fill = 0,
      names_glue = "{step}_{.value}"
    )
    df <- dplyr::filter(
      df, if_all(all_of(c("Before_count", "After_count")), \(y) y >= min_count)
    )
    df$count_change <- df$After_count - df$Before_count
    df$prop_change <- df$count_change / df$Before_count
    df
  }
)
#' @rdname summariseOverrep-methods
#' @importFrom dplyr bind_rows across summarise
#' @importFrom tidyselect all_of
#' @export
setMethod(
  "summariseOverrep", signature = "FastpDataList",
  function(
    x, min_count = 0, step = c("Before", "After"), vals = c("count", "rate"),
    fn = c("mean", "sum", "max"), by = c("reads", "sequence"), ...
  ){
    step <- match.arg(step, several.ok = TRUE)
    vals <- match.arg(vals, several.ok = TRUE)
    fun_list <- lapply(fn, match.fun)
    names(fun_list) <- fn
    df_list <- lapply(x, summariseOverrep, min_count = min_count)
    df <- bind_rows(df_list)
    cols <- apply(expand.grid(step, vals), 1, paste, collapse = "_")
    df <- summarise(
      df,
      n_samples = dplyr::n(),
      dplyr::across(all_of(cols), .fns = fun_list),
      .by = all_of(by)
    )
    df
  }
)
#' @importFrom dplyr across summarise
#' @importFrom tidyselect all_of
#' @rdname summariseOverrep-methods
#' @export
setMethod(
  "summariseOverrep", signature = "FastqcDataList",
  function(
    x, min_count = 0, vals = c("Count", "Percentage"),
    fn = c("mean", "sum", "max"), pattern = ".*", ...
  ){
    vals <- match.arg(vals, several.ok = TRUE)
    fun_list <- lapply(fn, match.fun)
    names(fun_list) <- fn
    df <- getModule(x, "Overrepresented_sequences")
    df <- df[grepl(pattern, df[["Possible_Source"]]),]
    df <- df[df[["Count"]] > min_count,]
    df <- summarise(
      df,
      n_samples = dplyr::n(),
      dplyr::across(all_of(vals), .fns = fun_list),
      .by = all_of(c("Sequence", "Possible_Source"))
    )
    df
  }
)
#' @importFrom dplyr across summarise
#' @importFrom tidyselect all_of
#' @rdname summariseOverrep-methods
#' @export
setMethod(
  "summariseOverrep", signature = "FastqcData",
  function(
    x, min_count = 0, vals = c("Count", "Percentage"),
    fn = c("mean", "sum", "max"), pattern = ".*", by = "Filename", ...
  ){
    vals <- match.arg(vals, several.ok = TRUE)
    fun_list <- lapply(fn, match.fun)
    names(fun_list) <- fn
    df <- getModule(x, "Overrepresented_sequences")
    df <- df[grepl(pattern, df[["Possible_Source"]]),]
    df <- df[df[["Count"]] > min_count,]
    if (nrow(df) == 0) return(df)
    summarise(
      df,
      n_sequences = dplyr::n(),
      dplyr::across(all_of(vals), .fns = fun_list),
      .by = all_of(by)
    )
  }
)
