suppressPackageStartupMessages({ library(tidyverse); library(readxl); library(decontam) })

org <- "Genus"

# get genus sheet from supplementary tables 1 and 2
# sheet 1 is raw phyla counts and  sheet 2 is raw genera counts
sheet_i <- 2

tax_df <- as.data.frame(read_xlsx("SUPPLEMENTARY_TABLE_S1.xlsx", sheet = sheet_i))
ntc_df <- as.data.frame(read_xlsx("SUPPLEMENTARY_TABLE_S2.xlsx", sheet = sheet_i))
metadata <- read_tsv("metadata.tsv", show_col_types = FALSE)

taxn_df <- bind_rows(tax_df, ntc_df)
taxn_df[is.na(taxn_df)] <- 0

indices <- match(taxn_df$Sample_Name, metadata$Samples)
stopifnot(!anyNA(indices))

conc_v <- as.numeric(metadata$Concentration[indices])
ntc_v  <- as.logical(metadata$NTC[indices])
stopifnot(!anyNA(conc_v), !anyNA(ntc_v), sum(ntc_v) == 1)

rownames(taxn_df) <- taxn_df$Sample_Name
taxn_mx <- as.matrix(taxn_df[, -1])

dir.create(org, showWarnings = FALSE)

# combined (primary method)
run_combined <- function(thr) {
  cont_df <- isContaminant(seqtab = taxn_mx, conc = conc_v, neg = ntc_v,
                           method = "combined", threshold = thr)
  list(cont = cont_df, keep = rownames(cont_df)[!(cont_df$contaminant %in% TRUE)])
}

# frequency only (NTC excluded)
smp_mx   <- taxn_mx[!ntc_v, , drop = FALSE]
smp_conc <- conc_v[!ntc_v]
run_freq <- function(thr) {
  cont_df <- isContaminant(seqtab = smp_mx, conc = smp_conc,
                           method = "frequency", threshold = thr)
  list(cont = cont_df, keep = rownames(cont_df)[!(cont_df$contaminant %in% TRUE)])
}

res_comb <- lapply(c(0, 0.10, 0.25, 0.50), run_combined)
names(res_comb) <- paste0("thr_", c(0, 0.10, 0.25, 0.50))

res_freq <- lapply(c(0.10, 0.25, 0.50), run_freq)
names(res_freq) <- paste0("freq_", c(0.10, 0.25, 0.50))

# strict (rule-based, no decontam)
ntc_taxa    <- colnames(taxn_mx)[colSums(taxn_mx[ntc_v, , drop = FALSE]) > 0]
strict_keep <- setdiff(colnames(taxn_mx), ntc_taxa)

res <- c(res_comb, res_freq, list(strict_NTC = list(cont = NULL, keep = strict_keep)))

# write per-scenario tables
for (nm in names(res)) {
  dir.create(file.path(org, nm), showWarnings = FALSE)
  if (!is.null(res[[nm]]$cont))
    write_tsv(rownames_to_column(res[[nm]]$cont, "taxon"),
              file.path(org, nm, "decontam_pvalues.tsv"))

  clean_df <- taxn_df[, c("Sample_Name", res[[nm]]$keep)]
  write_tsv(clean_df, file.path(org, nm, "taxa_wNTC.tsv"))     # NTC row kept, zero columns kept

  smp_clean <- clean_df[!ntc_v, ]
  smp_clean <- smp_clean[, c(TRUE, colSums(smp_clean[, -1, drop = FALSE]) > 0)]  # drop genera absent from all samples
  write_tsv(smp_clean, file.path(org, nm, "taxa_woNTC.tsv"))
}

# decontam summary
summary_df <- map_dfr(names(res), function(nm) {
  removed <- setdiff(colnames(taxn_mx), res[[nm]]$keep)
  p <- rowSums(smp_mx[, removed, drop = FALSE]) / rowSums(smp_mx) * 100
  tibble(
    scenario = nm,
    method = case_when(grepl("^freq", nm) ~ "frequency only",
                       nm == "strict_NTC" ~ "rule-based (no decontam)",
                       TRUE ~ "combined"),
    threshold = suppressWarnings(as.numeric(sub(".*_", "", nm))),
    genera_removed = length(removed),
    genera_removed_in_NTC = sum(removed %in% ntc_taxa),
    pct_reads_removed_min = min(p),
    pct_reads_removed_median = median(p),
    pct_reads_removed_max = max(p),
    NTC_reads_retained = sum(taxn_mx[ntc_v, res[[nm]]$keep, drop = FALSE]),
    removed_genera = paste(removed, collapse = "; "))
})
summary_df$threshold[summary_df$scenario == "strict_NTC"] <- NA
write_tsv(summary_df, file.path(org, "decontam_sensitivity_summary.tsv"))

# read depth (samples and NTC)
write_tsv(tibble(Sample_Name = rownames(taxn_mx), total_reads = rowSums(taxn_mx), is_NTC = ntc_v),
          file.path(org, "read_depth.tsv"))

# per-sample reads removed / retained, per scenario
per_sample <- map_dfr(names(res), function(nm) {
  removed <- setdiff(colnames(taxn_mx), res[[nm]]$keep)
  tibble(scenario = nm, Sample_Name = rownames(smp_mx),
         total_reads    = rowSums(smp_mx),
         retained_reads = rowSums(smp_mx[, res[[nm]]$keep, drop = FALSE]),
         frac_retained  = retained_reads / total_reads,
         pct_removed    = rowSums(smp_mx[, removed, drop = FALSE]) / total_reads * 100)
})
write_tsv(per_sample, file.path(org, "per_sample_pct_removed.tsv"))
per_sample |> group_by(scenario) |>
  summarise(min = min(frac_retained), median = median(frac_retained), max = max(frac_retained))

writeLines(capture.output(sessionInfo()), file.path(org, "sessionInfo.txt"))

pv_comb <- res$thr_0.1$cont |> rownames_to_column("taxon") |>
  select(taxon, p.freq, p.prev, p_combined = p)
pv_freq <- res$freq_0.1$cont |> rownames_to_column("taxon") |>
  select(taxon, p_freq_only = p)

pv_all <- pv_comb |> left_join(pv_freq, by = "taxon") |>
  mutate(in_NTC     = taxon %in% ntc_taxa,
         NTC_reads  = taxn_mx[ntc_v, taxon],
         in_samples = colSums(smp_mx)[taxon] > 0,
         untestable_combined = is.na(p_combined))
write_tsv(pv_all, file.path(org, "decontam_pvalues_all.tsv"))

# get the numbers used in the response
table(flagged = pv_all$p_combined < 0.1, in_NTC = pv_all$in_NTC, useNA = "ifany")
pv_all |> filter(in_NTC, is.na(p_combined)) |> pull(taxon) # the untestable NTC genera

summary(conc_v[!ntc_v]); conc_v[ntc_v] # what concentration did the NTC get?
arrange(tibble(Sample_Name = rownames(smp_mx), conc = conc_v[!ntc_v],
               reads = rowSums(smp_mx)), desc(conc)) # are the 3 high samples also the deep ones?
cor(log10(conc_v[!ntc_v]), log10(rowSums(smp_mx)), method = "spearman")
