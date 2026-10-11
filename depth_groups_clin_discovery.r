source("common.r") # loads libraries, org, scenarios, load_tss and dat

# depth groups (based on unfiltered reads)
tot0 <- setNames(dat$thr_0$totals, dat$thr_0$counts$Sample_Name)
is_ntc0 <- grepl("NTC", names(tot0), ignore.case = TRUE)
ntc_reads <- tot0[is_ntc0]                       # ~747

depth_df <- tibble(Sample_Name = names(tot0)[!is_ntc0],
                   reads = unname(tot0[!is_ntc0])) |>
  mutate(group = if_else(reads < ntc_reads, "low (< NTC)", "adequate"),
         ntc_to_sample = ntc_reads / reads)    # worst-case background fraction
write_tsv(depth_df, file.path(org, "depth_groups.tsv"))
count(depth_df, group)

# dominant genus per scenario vs unfiltered
dom <- map_dfr(scenarios, function(s) {
  r <- dat[[s]]$rel |> filter(!grepl("NTC", Sample_Name, ignore.case = TRUE))
  m <- as.matrix(r[, -1])
  tibble(scenario = s, Sample_Name = r$Sample_Name,
         top_genus = colnames(m)[max.col(m, ties.method = "first")],
         top_ra = apply(m, 1, max))
})
base_top <- dom |> filter(scenario == "thr_0") |> select(Sample_Name, base_genus = top_genus)

dom_cmp <- dom |> left_join(base_top, by = "Sample_Name") |>
  left_join(depth_df, by = "Sample_Name") |>
  mutate(same = top_genus == base_genus)

ret <- map_dfr(scenarios, function(s)
  tibble(scenario = s, Sample_Name = dat[[s]]$counts$Sample_Name, retained_reads = dat[[s]]$totals)) |>
  filter(!grepl("NTC", Sample_Name, ignore.case = TRUE))

dom_cmp <- dom_cmp |>
  left_join(ret, by = c("scenario", "Sample_Name")) |>
  mutate(interpretable = retained_reads >= 100) # percentages of tiny totals are not meaningful

dom_cmp |> filter(scenario == "strict_NTC") |>
  select(Sample_Name, reads, retained_reads, base_genus, top_genus, top_ra, same, interpretable) |>
  print(n = Inf)

dom_cmp |> filter(scenario != "thr_0") |>
  group_by(scenario, group) |>
  summarise(n = n(), n_same = sum(same),
            n_interpretable = sum(interpretable),
            n_same_interpretable = sum(same & interpretable), .groups = "drop") |>
  print(n = Inf)

# samples whose dominant genus changed (any scenario)
dom_cmp |> filter(!same) |>
  select(scenario, Sample_Name, reads, retained_reads, base_genus, top_genus, interpretable) |>
  print(n = Inf)

# Bray-Curtis vs unfiltered, by depth group
base <- dat$thr_0$rel |> filter(!grepl("NTC", Sample_Name, ignore.case = TRUE))
bc <- map_dfr(setdiff(scenarios, "thr_0"), function(s) {
  r <- as.data.frame(dat[[s]]$rel)
  r <- r[match(base$Sample_Name, r$Sample_Name), , drop = FALSE]
  miss <- setdiff(colnames(base)[-1], colnames(r)); r[miss] <- 0
  r <- r[, colnames(base)]
  tibble(scenario = s, Sample_Name = base$Sample_Name,
         bray_curtis = 1 - rowSums(pmin(as.matrix(base[, -1]), as.matrix(r[, -1]))))
}) |> left_join(depth_df, by = "Sample_Name")

bc |> group_by(scenario, group) |>
  summarise(n = n(), median = median(bray_curtis), max = max(bray_curtis), .groups = "drop") |>
  print(n = Inf)

# sample-clinically detected-expected bacteria
chk <- tribble(~Sample_Name, ~genus,
  "BOS09","Klebsiella", "BOS11","Klebsiella", "BOS16","Enterococcus",
  "BOS21","Mycobacterium", "BOS23","Neisseria")

cnt <- dat$thr_0$counts

clin <- chk |> rowwise() |> mutate(
  reads_in_sample = if (genus %in% colnames(cnt)) cnt[[genus]][cnt$Sample_Name == Sample_Name] else NA,
  sample_total    = cnt[rowSums(cnt[,-1]) > 0 & cnt$Sample_Name == Sample_Name, ] |> select(-1) |> sum(),
  reads_in_NTC    = if (genus %in% colnames(cnt)) cnt[[genus]][grepl("NTC", cnt$Sample_Name)] else NA)

write_tsv(clin, file.path(org, "sequencing_discovery.tsv"))

keep_tbl <- ret |> left_join(select(depth_df, Sample_Name, reads), by = "Sample_Name") |>
  mutate(frac_retained = retained_reads / reads)
write_tsv(keep_tbl, file.path(org, "fraction_reads_retained.tsv"))

keep_tbl |> group_by(scenario) |>
  summarise(min = min(frac_retained), median = median(frac_retained), max = max(frac_retained))

set.seed(1)
cnt0 <- dat$thr_0$counts |> filter(!grepl("NTC", Sample_Name, ignore.case = TRUE))
m0 <- as.matrix(cnt0[, -1]); rownames(m0) <- cnt0$Sample_Name
stab <- map_dfr(rownames(m0), function(sn) {
  x <- m0[sn, ]; top <- which.max(x)
  draws <- rmultinom(1000, size = sum(x), prob = x / sum(x))
  tibble(Sample_Name = sn, reads = sum(x), top_genus = names(x)[top],
         top_ra = max(x) / sum(x),
         p_top_stays_top = mean(apply(draws, 2, which.max) == top))
}) |> left_join(select(depth_df, Sample_Name, group), by = "Sample_Name") |> arrange(reads)
print(stab, n = Inf)
write_tsv(stab, file.path(org, "dominant_genus_stability.tsv"))

wide_dom <- dom_cmp |> select(scenario, Sample_Name, top_genus) |>
  pivot_wider(names_from = scenario, values_from = top_genus, names_prefix = "top_genus_")
wide_ret <- ret |>
  pivot_wider(names_from = scenario, values_from = retained_reads, names_prefix = "retained_reads_")
wide_bc  <- bc |> select(scenario, Sample_Name, bray_curtis) |>
  pivot_wider(names_from = scenario, values_from = bray_curtis, names_prefix = "bray_curtis_")

s6 <- depth_df |>
  left_join(select(stab, Sample_Name, top_ra, p_top_stays_top), by = "Sample_Name") |>
  left_join(wide_dom, by = "Sample_Name") |>
  left_join(wide_ret, by = "Sample_Name") |>
  left_join(wide_bc,  by = "Sample_Name")
write_tsv(s6, file.path(org, "S6_per_sample.tsv"))
