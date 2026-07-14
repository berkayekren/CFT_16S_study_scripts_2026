# CFT analysis R script #2

suppressPackageStartupMessages({
  library(tidyverse)
})

# import cleaned taxa for relative abundance calculations (TSS)
org <- "Phylum" # Genus
taxn_clean_df <- read_tsv(file.path(org, "taxa_wNTC.tsv"))

# tss calculation
row_totals <- rowSums(taxn_clean_df[, -1])
tss_df <- taxn_clean_df
tss_df[, -1] <- taxn_clean_df[, -1] / row_totals
tss_df[is.na(tss_df)] <- 0

# export the tss file
write_tsv(tss_df, file.path(org, "tss_wNTC.tsv"))
