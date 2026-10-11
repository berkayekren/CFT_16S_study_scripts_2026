suppressPackageStartupMessages({ library(tidyverse); library(RColorBrewer) })

org <- "Genus"
scenarios <- c("thr_0", "thr_0.1", "thr_0.25", "thr_0.5",
               "freq_0.1", "freq_0.25", "freq_0.5", "strict_NTC")

# load counts and compute relative abundance (after removal)
load_tss <- function(s) {
  cnt <- read_tsv(file.path(org, s, "taxa_wNTC.tsv"), show_col_types = FALSE)
  tot <- rowSums(cnt[, -1])
  is_ntc <- grepl("NTC", cnt$Sample_Name, ignore.case = TRUE)
  stopifnot(all(tot[!is_ntc] > 0))        # real samples must keep some reads
  rel <- cnt
  rel[, -1] <- cnt[, -1] / pmax(tot, 1)   # NTC row stays 0 instead of NaN
  list(counts = cnt, rel = rel, totals = tot)
}
dat <- map(set_names(scenarios), load_tss)
