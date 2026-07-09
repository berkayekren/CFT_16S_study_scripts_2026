# CFT analysis R script #ST

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

  my_colours <- c("#666666", "#ffa500", "#008b00", "#0000cd", "#570826", brewer.pal(7, "Dark2")[c(3, 1, 2, 4, 5, 6, 7)], "#773935", "#900C3F", brewer.pal(12, "Set3")[-c(2, 9)], brewer.pal(8, "Set1")[1:8], brewer.pal(7, "Accent")[1:7], brewer.pal(10, "Paired"))

  org_col_name <- paste0("Organisms (", my_tax, ifelse(seq_type == "", "", "-"), seq_type, "):")

  tax_st <- rat_ss %>%
    rename(Samples = 1, Taxon = 2, RA = 3, Group = 4) %>%
    mutate(
      Samples = factor(Samples, levels = unique(Samples), ordered = TRUE),
      RA = as.numeric(RA)
    )

  if (my_tax != "Phylum") {
    tax_st <- tax_st %>%
      mutate(
        # Strip " spp." for the lookup
        Clean_Tax = str_remove(Taxon, " spp."),
        # Instantly look up the corresponding Phylum from otu_total
        Phylum_Match = otu_total$Phylum[match(Clean_Tax, otu_total[[my_tax]])],
        # Prepend the Phylum if it's not "Others" and a match was found
        Taxon = ifelse(Taxon != " Others" & !is.na(Phylum_Match),
                       paste0(Phylum_Match, "_", Taxon),
                       Taxon)
      ) %>%
      select(-Clean_Tax, -Phylum_Match) # drop the temporary columns
  }

  tax_st <- tax_st %>% rename(!!sym(org_col_name) := Taxon)
  n_taxa <- n_distinct(tax_st[[org_col_name]])

  if (n_taxa <= 40) {
    col_num <- ifelse(n_taxa <= 25, 1, 2)
    s_cats  <- ifelse(n_taxa <= 25, 0, 4)
    plot_st <- ggplot(tax_st, aes(x = Samples, y = RA, fill = .data[[org_col_name]])) +
      geom_bar(position = "stack", stat = "identity", width = 0.95) +
      scale_fill_manual(values = my_colours, name = org_col_name) +
      guides(fill = guide_legend(ncol = col_num)) +
      labs(
        y = "Relative Abundance (%)",
        x = "",
        title = paste0(tax_name, ifelse(seq_type == "", "", paste0(" (", seq_type, ")")), "\n")
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    save_path <- paste0("StackedAbundance_", my_tax, ".tiff")
    ggsave(filename = save_path, limitsize = FALSE, plot = plot_st, units = "in", width = (round(s_all / 6, 0) + 7 + s_cats), height = 8, dpi = 600, compression = "lzw")}
