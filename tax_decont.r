# CFT analysis R script #1

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(decontam)
})


# import & create matrices and lists
meta_tax <- read_tsv(file.path("meta_tax.tsv"))
metadata <- read_tsv(file.path("metadata.tsv"))
conc_v <- unlist(metadata$concentrations)

# get sample specifications
neg_v <- unlist(metadata$NTC)
rownames(metadata) <- metadata$Samples
meta_tax_df <- as.data.frame(meta_tax)
rownames(meta_tax_df) <- paste0("TAX", seq_len(nrow(meta_tax_df)))

# transpose matrix for decontam
meta_mx <- as.matrix(t(meta_tax_df[, -1]))

# find possible contaminants and remove them from the main meta_tax
cont_df <- isContaminant(seqtab = meta_mx, conc = conc_v, neg = neg_v, method = "combined")

# original organism names for later reference
cont_df$Orgs <- meta_tax_df[, 1]

# filter the contaminants
cont_tbl <- cont_df[which(cont_df$contaminant), ]
meta_tax_clean <- meta_tax_df[!(rownames(meta_tax_df) %in% rownames(cont_tbl)), ]

# remove NTC from the analysis list
ntc_cols <- grep("NTC", colnames(meta_tax_clean))
if (length(ntc_cols) > 0) {
  meta_tax_clean <- meta_tax_clean[, -ntc_cols]
}

# keep the row names and export
meta_tax_clean %>%
  rownames_to_column(var = "TAX_ID") %>%
  write_tsv("meta_tax_decon.tsv")
