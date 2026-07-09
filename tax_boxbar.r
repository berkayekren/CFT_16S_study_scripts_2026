# CFT analysis R script #BBP

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
})

# import necessary files
meta_tax <- read_tsv(file.path("meta_tax_decon.tsv"))
metadata <- read_tsv(file.path("metadata.tsv"))
otu_total <- read_excel("~/all_taxa_file.xlsx", sheet = 1) %>%
  distinct() %>%
  arrange(pick(2))

# get all samples without ntc
s_all <- nrow(metadata) - 1
tax_name <- "Phylum" # or "Genus"

# pivot the wide table to a long format
long_data <- rel_abund %>%
  rename(Taxon = 1) %>%
  pivot_longer(
    cols = -Taxon,
    names_to = "V1",
    values_to = "V3"
  ) %>%
  rename(V2 = Taxon)

# calculate the global average ("BOS-AVE")
ave_data <- long_data %>%
  group_by(V2) %>%
  summarize(V3 = mean(V3, na.rm = TRUE), .groups = "drop") %>%
  mutate(V1 = "BOS-AVE")

# combine, apply the 1% limit, and aggregate
otu_stack <- bind_rows(long_data, ave_data) %>%
  mutate(V2 = ifelse(V3 < 1, " Others", V2)) %>%
  group_by(V1, V2) %>%
  summarize(V3 = sum(V3, na.rm = TRUE), .groups = "drop") %>%
  arrange(V1, V2)

# prepare for plots
my_colours <- c("#666666", "#ffa500", "#008b00", "#0000cd", "#570826", brewer.pal(7, "Dark2")[c(3, 1, 2, 4, 5, 6, 7)], "#773935", "#900C3F", brewer.pal(12, "Set3")[-c(2, 9)], brewer.pal(8, "Set1")[1:8], brewer.pal(7, "Accent")[1:7], brewer.pal(10, "Paired"))
tax_la <- otu_stack[otu_stack[, 3] > 1 && ! " Others " %in% otu_stack[, 2], ]
tax_la <- tax_la[order(tax_la[, 3], decreasing = T), ]

# log plot variance
tax_ls <- cbind(tax_la[, 1:2], log10(tax_la[, 3]))
tax_ls[, 3] <- as.numeric(tax_ls[, 3])
colnames(tax_ls) <- c("Samples", "Organisms", "log10")
plot_ls <- ggplot(data = tax_ls, aes(x = Organisms, y = log10)) +
  stat_boxplot(geom = "errorbar", linewidth = 0.4, coef = 1.5) +
  scale_fill_manual(values = my_colours[-1]) +
  geom_boxplot(outlier.colour = NA, aes(fill = Organisms), colour = "black", linewidth = 0.2) +
  labs(title = paste0(tax_name, ifelse(seq_type == "", "", paste0(" (", seq_type, ")")), "\n"), y = "Relative Abundance (log10)", x = "") +
  theme(panel.grid.major = element_line(colour = "black", linewidth = 0.1), panel.background = element_rect(fill = "white"), legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), plot.title = element_text(size = 10), axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.9)) +
  coord_flip(clip = "off")
ggsave(filename = paste0("RelLogVertBox_", my_tax, "_all.tiff"), plot = plot_la, units = "in", width = 8, height = 7, dpi = 600, compression = "lzw")
