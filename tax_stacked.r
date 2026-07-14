# CFT analysis R script #3

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
})

# load the data
org <- "Phylum" # "Genus"
clean_counts <- read_tsv(file.path(org, "taxa_wNTC.tsv"))
tss_wide <- read_tsv(file.path(org, "tss_wNTC.tsv"))

# create the long-table from wide
tss_long <- pivot_longer(
  data = tss_wide,
  cols = -Sample_Name,
  names_to = "Organism",
  values_to = "Relative_Abundance"
)

# filter the <0.01
tss_long$Organism[tss_long$Relative_Abundance < 0.01] <- "Others"
plot_data_collapsed <- tss_long %>%
  group_by(Sample_Name, Organism) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance), .groups = "drop")

# create GAP for plot design
sample_totals <- rowSums(clean_counts[, -1])
label_map <- data.frame(
  Sample_Name = clean_counts$Sample_Name,
  X_Label = paste0(clean_counts$Sample_Name, "\n(n=", sample_totals, ")"),
  stringsAsFactors = FALSE
)

# get raw read counts
is_ntc <- grepl("NTC", label_map$Sample_Name, ignore.case = TRUE)
smp_labels <- label_map$X_Label[!is_ntc]
ntc_labels <- label_map$X_Label[is_ntc]

# set up the plot design
gap_size <- 1.4
smp_positions <- seq_along(smp_labels)
ntc_positions <- (length(smp_labels) + gap_size) + (seq_along(ntc_labels) - 1)
position_df <- data.frame(
  X_Label = c(smp_labels, ntc_labels),
  X_Numeric = c(smp_positions, ntc_positions),
  stringsAsFactors = FALSE
)

# set the sample order in plot
plot_data <- left_join(plot_data_collapsed, label_map, by = "Sample_Name")
plot_data <- left_join(plot_data, position_df, by = "X_Label")

# set the "Other" section at the top
clean_orgs_only <- sort(setdiff(unique(plot_data$Organism), "Others"))
org_levels <- c("Others", clean_orgs_only)
plot_data$Organism <- factor(plot_data$Organism, levels = org_levels)

# set the colors
my_colours <- c(
  "#E0E0E0", "#ffa500", "#008b00", "#0000cd", "#570826",
  brewer.pal(7, "Dark2")[c(3, 1, 2, 4, 5, 6, 7)],
  "#773935", "#900C3F",
  brewer.pal(12, "Set3")[-c(2, 9)],
  brewer.pal(8, "Set1")[1:8],
  brewer.pal(7, "Accent")[1:7],
  brewer.pal(10, "Paired")
)

# run plot
st_plot <- ggplot(plot_data, aes(x = X_Numeric, y = Relative_Abundance, fill = Organism)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.2) +
  scale_x_continuous(breaks = position_df$X_Numeric, labels = position_df$X_Label) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = my_colours) +
  labs(
    x = "Sample ID (Clean Post-Decontam Reads)",
    y = "Relative Abundance (%)",
    fill = paste0("Organisms: (", org, ")")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.text.x = element_text(face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

# save the plot
ggsave(
  filename = paste0(org, "/Taxonomic_Composition_", org, ".png"),
  plot = st_plot,
  width = 1.5 + nrow(clean_counts) + ceiling(length(org_levels) / 16),
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)
