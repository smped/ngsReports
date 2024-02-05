data <- "~/TKI/DBPAInT/output/fastp/23-00108_S1_L01_fastp.json" %>%
  jsonlite::fromJSON()

data <- "inst/extdata/fastp.json.gz" %>%
  jsonlite::fromJSON()

c(
  "read1_before_filtering", "read2_before_filtering",
  "read1_after_filtering", "read2_after_filtering"
  ) %>%
  sapply(
    \(x) unlist(data[[x]]$overrepresented_sequences),
    simplify = FALSE
  ) %>%
  lapply(as_tibble, rownames = "sequence") %>%
  bind_rows(.id = "step") %>%
  separate(step, c("reads", "filtering"), extra = "merge") %>%
  pivot_wider(names_from = "filtering", values_from = "value", values_fill = 0) %>%
  arrange(desc(after_filtering))


%>%
  .[["read1_before_filtering"]] %>%
  .[["overrepresented_sequences"]] %>%
  unlist() %>%
  as_tibble(rownames = "sequence") %>%
  dplyr::rename(count = value) %>%
  arrange(desc(count))
