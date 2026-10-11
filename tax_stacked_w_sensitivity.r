source("common.r") # loads libraries, org, scenarios, load_tss and dat

# one fixed palette across all scenarios
my_colours <- c(
  "#ffa500", "#008b00", "#0000cd", "#570826",
  brewer.pal(7, "Dark2")[c(3, 1, 2, 4, 5, 6, 7)],
  "#773935", "#900C3F",
  brewer.pal(12, "Set3")[-c(2, 9)],
  brewer.pal(8, "Set1")[1:8],
  brewer.pal(7, "Accent")[1:7],
  brewer.pal(10, "Paired")
)

all_long <- map_dfr(dat, function(d)
  d$rel |> filter(!grepl("NTC", Sample_Name)) |>
    pivot_longer(-Sample_Name, names_to = "Organism", values_to = "ra"),
  .id = "scenario")

mean_order <- all_long |> group_by(Organism) |> summarise(m = mean(ra)) |>
  arrange(desc(m)) |> pull(Organism)

dominant <- all_long |> group_by(scenario, Sample_Name) |>
  slice_max(ra, n = 1, with_ties = FALSE) |> pull(Organism)

shown <- intersect(mean_order, unique(c(mean_order[1:35], dominant)))
message(length(shown), " genera coloured")
stopifnot(length(shown) <= length(my_colours))

pal <- c(Others = "#E0E0E0", setNames(my_colours[seq_along(shown)], shown))
stopifnot(!anyNA(pal))

# plotting function
ntc_ref <- unname(dat$thr_0$totals[grepl("NTC", dat$thr_0$counts$Sample_Name, ignore.case = TRUE)])
stopifnot(length(ntc_ref) == 1)

plot_scenario <- function(s) {
  d <- dat[[s]]
  long <- d$rel |>
    pivot_longer(-Sample_Name, names_to = "Organism", values_to = "Relative_Abundance") |>
    mutate(Organism = if_else(Organism %in% shown, Organism, "Others")) |>
    group_by(Sample_Name, Organism) |>
    summarise(Relative_Abundance = sum(Relative_Abundance), .groups = "drop")

  label_map <- tibble(Sample_Name = d$counts$Sample_Name,
                      reads = d$totals,
                      X_Label = paste0(Sample_Name, "\n(n=", d$totals, ")"))
  is_ntc <- grepl("NTC", label_map$Sample_Name, ignore.case = TRUE)
  n_smp <- sum(!is_ntc)
  pos <- tibble(X_Label = c(label_map$X_Label[!is_ntc], label_map$X_Label[is_ntc]),
                X_Numeric = c(seq_len(n_smp), n_smp + 1.4 + seq_len(sum(is_ntc)) - 1))

  low <- label_map$reads < ntc_ref & !is_ntc
  lab_cols <- setNames(ifelse(low, "red3", "black"), label_map$X_Label)[pos$X_Label]

  pd <- long |> left_join(label_map, by = "Sample_Name") |> left_join(pos, by = "X_Label") |>
    mutate(Organism = factor(Organism, levels = c("Others", shown)))
  stopifnot(!anyNA(pd$Organism), !anyNA(pd$X_Numeric))

  ggplot(pd, aes(X_Numeric, Relative_Abundance, fill = Organism)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.2) +
    scale_x_continuous(breaks = pos$X_Numeric, labels = pos$X_Label) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = pal, breaks = names(pal), limits = names(pal),
                      drop = FALSE, na.value = "black") +
    labs(title = s, x = "Sample ID (retained reads after this scenario)",
         y = "Relative Abundance (%)", fill = paste0("Organisms: (", org, ")")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
          axis.text.x = element_text(face = "bold", size = 8, color = unname(lab_cols)),
          axis.title = element_text(face = "bold", size = 11),
          legend.title = element_text(face = "bold", size = 10))
}

walk(scenarios, function(s) {
  ggsave(file.path(org, s, "Taxonomic_Composition.png"), plot_scenario(s),
         width = 1.5 + nrow(dat[[s]]$counts) + 2, height = 7.5, dpi = 600, bg = "white")
  message("saved plot for ", s)
})

# combined supplementary figure: vector PDF, one shared legend, panel letters A-H
library(patchwork)

plots <- map(scenarios, plot_scenario)
combined <- wrap_plots(plots, ncol = 1) + plot_annotation(tag_levels = "A")

ggsave(file.path(org, "Supplementary_Figure_S1.pdf"), combined,
       width = 32, height = 7.5 * length(plots), units = "in",
       limitsize = FALSE, device = cairo_pdf)
