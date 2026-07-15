# CFT analysis R script #1

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(decontam)
})

# import & create matrices and lists
org <- "Phylum" # "Phylum"
tax_df <- as.data.frame(read_xlsx(file.path("SUPPLEMENTARY_TABLE1.xlsx"), sheet=1)) # sheet=1 for phylum; 2 for genus
ntc_df <- as.data.frame(read_xlsx(file.path("SUPPLEMENTARY_TABLE2.xlsx"), sheet=1)) # sheet=1 for phylum; 2 for genus
metadata <- read_tsv(file.path("metadata.tsv"))

# get sample specifications
taxn_df <- bind_rows(tax_df , ntc_df)
taxn_df[is.na(taxn_df)] <- 0
indices <-match(taxn_df$Sample_Name, metadata$Samples)
conc_v <- unlist(metadata$Concentration[indices])
ntc_v <- unlist(metadata$NTC[indices])
rownames(taxn_df) <- taxn_df$Sample_Name
taxn_mx <- as.matrix(taxn_df[, -1])

# find possible contaminants and remove them from the main taxn_df
cont_df <- isContaminant(seqtab = taxn_mx, conc = conc_v, neg = ntc_v, method = "combined")
clean_taxa_names <- rownames(cont_df)[cont_df$contaminant == FALSE]
columns_to_keep <- c("Sample_Name", clean_taxa_names)
taxn_clean_df <- taxn_df[, columns_to_keep]
tax_clean_df <- taxn_clean_df[! taxn_clean_df$Sample_Name %in% "NTC"]

# create the directories
dir.create(org)

# export cleaned tax with and without NTC
write_tsv(taxn_clean_df, file.path(org, "taxa_wNTC.tsv"))
write_tsv(tax_clean_df, file.path(org, "taxa_woNTC.tsv"))
