# CFT 16S Study Scripts (2026)

Analysis scripts for the article **"Microbiota Profiling of Cerebrospinal Fluid in Bacterial Meningitis Using 16S rRNA Gene Sequencing"** (16S rRNA gene sequencing, ONT platform).

The scripts identify putative contaminants with [decontam](https://github.com/benjjneb/decontam) and test how sensitive the genus-level composition is to the filtering choice and to sequencing depth.


## NCBI Uploads & IDs

16 of 26 samples and NTC were uploaded to the NCBI SRA Database. The given BioProject ID is `PRJNA1491059`. The SRA IDs assigned are given below.

| ID | SRA Number |
| BOS12	| SRR39498336 |
| BOS13	| SRR39498337 |
| BOS14	| SRR39498338 |
| BOS15	| SRR39498339 |
| BOS16	| SRR39498340 |
| BOS17	| SRR39498341 |
| BOS18	| SRR39498328 |
| BOS19	| SRR39498329 |
| BOS20	| SRR39498330 |
| BOS21	| SRR39498331 |
| BOS22	| SRR39498332 |
| BOS23	| SRR39498333 |
| BOS24	| SRR39498334 |
| BOS25	| SRR39498335 |
| BOS26	| SRR39498342 |
| NTC	| SRR39498343 |

The samples BOS01 to 11 will be uploaded after recovery.


## Scenarios

All scenarios are run on genus-level raw read counts of the study samples plus one no-template control (NTC).

| Scenario | Method |
|---|---|
| `thr_0` | No filtering |
| `thr_0.1` | decontam, `combined` method, threshold 0.1 (primary, decontam default) |
| `thr_0.25`, `thr_0.5` | decontam, `combined` method, more permissive thresholds |
| `freq_0.1`, `freq_0.25`, `freq_0.5` | decontam, `frequency` method alone (study samples only, NTC not used) |
| `strict_NTC` | Rule-based: removes every genus detected in the NTC (no decontam) |

Relative abundances are calculated by total-sum scaling of the reads retained in each scenario.

## Scripts

Run from the repository root, in this order:

| Order | File | Purpose |
|---|---|---|
| 1 | `tax_decontam_w_sensitivity.r` | Runs decontam and the strict scenario; writes per-scenario count tables, p-values, read depth and the scenario summary |
| (sourced) | `common.r` | Shared setup: libraries, scenario names, loading of count tables and total-sum scaling |
| 2 | `tax_stacked_w_sensitivity.r` | Stacked bar charts per scenario and the combined Supplementary Figure S1 |
| 3 | `depth_groups_clin_discovery.r` | Depth groups (samples below or above the NTC read count), dominant-genus comparison across scenarios, Bray-Curtis, resampling stability, per-sample table |

```bash
Rscript tax_decontam_w_sensitivity.r
Rscript tax_stacked_w_sensitivity.r
Rscript depth_groups_clin_discovery.r
```

## Input files

- `SUPPLEMENTARY_TABLE_S1.xlsx`: raw read counts per sample (sheet 1 phylum, sheet 2 genus)
- `SUPPLEMENTARY_TABLE_S2.xlsx`: raw read counts of the NTC (same sheet layout)
- `metadata.tsv`: columns `Samples`, `Concentration` (pre-sequencing DNA concentration) and `NTC` (TRUE for the control)

decontam is run on the **genus** sheet only.

## Output (`./` - project root)
- `sessionInfo.txt`: R and package versions
- `Genus/<scenario>/taxa_wNTC.tsv`, `taxa_woNTC.tsv`: count tables after filtering, with and without the NTC row
- `Genus/<scenario>/decontam_pvalues.tsv`: decontam p-values (not written for `strict_NTC`)
- `Genus/<scenario>/Taxonomic_Composition.png`: stacked bar chart
- `Genus/decontam_sensitivity_summary.tsv`: genera and reads removed per scenario
- `Genus/decontam_pvalues_all.tsv`: combined and frequency-only p-values per genus, with an NTC flag
- `Genus/read_depth.tsv`, `per_sample_pct_removed.tsv`, `fraction_reads_retained.tsv`: read depth and reads retained per sample
- `Genus/depth_groups.tsv`, `dominant_genus_stability.tsv`, `S6_per_sample.tsv`: depth-group and stability results
- `Genus/sequencing_discovery.tsv`: read counts of selected genera per sample
- `Genus/Supplementary_Figure_S1.pdf`: all scenarios in one figure


## Requirements

R with `tidyverse`, `readxl`, `decontam` (v1.24.0), `RColorBrewer` and `patchwork`. See `Genus/sessionInfo.txt` for the exact versions used.

## Notes

- Only one NTC was available, and it was sequenced with a single sequencing batch. The prevalence test of decontam therefore has very limited resolution, which is why the additional scenarios are provided.
- Samples with fewer reads than the NTC (747 reads) are marked in red in the figures and should be interpreted as low-confidence.
- `isContaminant()` is deterministic. The only random step is the resampling in `depth_groups_clin_discovery.r` (`set.seed(1)`).

## Citation

Article: "Microbiota Profiling of Cerebrospinal Fluid in Bacterial Meningitis Using Full-Length 16S rRNA Gene Amplicon Sequencing"
