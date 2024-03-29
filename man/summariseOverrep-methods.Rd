% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summariseOverrep.R
\name{summariseOverrep}
\alias{summariseOverrep}
\alias{summariseOverrep,FastpData-method}
\alias{summariseOverrep,FastpDataList-method}
\alias{summariseOverrep,FastqcDataList-method}
\alias{summariseOverrep,FastqcData-method}
\title{Summarise Overrepresented Sequences}
\usage{
summariseOverrep(x, ...)

\S4method{summariseOverrep}{FastpData}(x, step = c("Before", "After"), min_count = 0, ...)

\S4method{summariseOverrep}{FastpDataList}(
  x,
  min_count = 0,
  step = c("Before", "After"),
  vals = c("count", "rate"),
  fn = c("mean", "sum", "max"),
  by = c("reads", "sequence"),
  ...
)

\S4method{summariseOverrep}{FastqcDataList}(
  x,
  min_count = 0,
  vals = c("Count", "Percentage"),
  fn = c("mean", "sum", "max"),
  pattern = ".*",
  ...
)

\S4method{summariseOverrep}{FastqcData}(
  x,
  min_count = 0,
  vals = c("Count", "Percentage"),
  fn = c("mean", "sum", "max"),
  pattern = ".*",
  by = "Filename",
  ...
)
}
\arguments{
\item{x}{An object of a suitable class}

\item{...}{Not used}

\item{step}{Can be 'Before', 'After' or both to obtain data from the
Before_filtering or After_filtering modules}

\item{min_count}{Filter sequences with counts less than this value, both
before and after filtering}

\item{vals}{Values to use for creating summaries across multiple files. For
FastpDataList objects these can be "count" and/or "rate", whilst for
FastqcDataList objects these values can be "Count" and/or "Percentage"}

\item{fn}{Functions to use when summarising values across multiple files}

\item{by}{character vector of columns to summarise by. See \link[dplyr:summarise]{dplyr::summarise}}

\item{pattern}{Regular expression to filter the Possible_Source column by}
}
\value{
A tibble

Tibble columns will vary between Fastp*, FastqcDataList and FastqcData
objects. Calling this function on list-type objects will attempt to
summarise the presence each over-represented sequence across all files.

In particular, FastqcData objects will provide the requested summary
statistics across all sequences within a file
}
\description{
Summarise the Overrepresented sequences found in one or more QC files
}
\details{
This function prepares a useful summary of all over-represented sequences
as reported by either fastp or FastQC
}
\examples{
## For operations on a FastpData object
f <- system.file("extdata/fastp.json.gz", package = "ngsReports")
fp <- FastpData(f)
summariseOverrep(fp, min_count = 100)

## Applying the function to a FastqcDataList
packageDir <- system.file("extdata", package = "ngsReports")
fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
fdl <- FastqcDataList(fl)
summariseOverrep(fdl)

# An alternative viewpoint can be obtained using
fdl |> lapply(summariseOverrep) |> dplyr::bind_rows()


}
