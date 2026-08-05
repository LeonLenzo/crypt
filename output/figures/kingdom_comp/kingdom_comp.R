#!/usr/bin/env Rscript
# figure/kingdom_comp/kingdom_comp.R — secondary kingdom composition by mode
# Run from crypt/: Rscript figure/kingdom_comp/kingdom_comp.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
})

RUNS_TSV <- "output/02_filter_runs/data/runs.tsv"
OUT       <- "output/figures/kingdom_comp/kingdom_comp"

KINGDOM_COLS <- c(
  Fungi     = "#d35400",
  Viruses   = "#c0392b",
  Bacteria  = "#2980b9",
  Oomycota  = "#8e44ad",
  Nematoda  = "#27ae60",
  Unknown   = "#95a5a6"
)

crypt <- read.delim(RUNS_TSV, stringsAsFactors = FALSE) |>
  filter(biosample_representative == "True")

# ── Panel A: kingdom breakdown of secondary detections per mode ────────────────
coinf <- crypt |>
  filter(co_infection_flag != "single", secondary_kingdoms != "") |>
  select(mode, secondary_kingdoms) |>
  mutate(secondary_kingdoms = str_remove_all(secondary_kingdoms,
                                             "\\s*\\([^)]+\\)")) |>
  separate_rows(secondary_kingdoms, sep = ";\\s*") |>
  filter(secondary_kingdoms != "", !is.na(secondary_kingdoms)) |>
  mutate(
    kingdom   = str_trim(secondary_kingdoms),
    mode_label = ifelse(mode == "mal", "MAL", "HAL")
  ) |>
  count(mode_label, kingdom) |>
  group_by(mode_label) |>
  mutate(
    prop  = n / sum(n),
    kingdom = factor(kingdom, levels = names(KINGDOM_COLS))
  ) |>
  ungroup()

cat("Secondary kingdom counts:\n")
print(as.data.frame(coinf |> select(mode_label, kingdom, n, prop)))

p_a <- ggplot(coinf,
              aes(x = mode_label, y = prop, fill = kingdom)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(prop >= 0.03,
                               paste0(round(prop * 100), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand  = c(0, 0)) +
  scale_fill_manual(values = KINGDOM_COLS, name = "Secondary kingdom",
                    drop = FALSE) +
  labs(x = NULL, y = "Proportion of secondary detections",
       title = "Secondary kingdom composition") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

# ── Panel B: co-infection flag breakdown per mode ─────────────────────────────
flag_dat <- crypt |>
  mutate(
    mode_label        = ifelse(mode == "mal", "MAL", "HAL"),
    co_infection_flag = factor(co_infection_flag,
                               levels = c("multi_kingdom", "multi_species", "single"))
  ) |>
  count(mode_label, co_infection_flag) |>
  group_by(mode_label) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

FLAG_COLS <- c(
  single        = "#dfe6e9",
  multi_species = "#e67e22",
  multi_kingdom = "#c0392b"
)
FLAG_LABS <- c(
  single        = "Single",
  multi_species = "Multi-species",
  multi_kingdom = "Multi-kingdom"
)

p_b <- ggplot(flag_dat,
              aes(x = mode_label, y = prop, fill = co_infection_flag)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(prop >= 0.01,
                               paste0(round(prop * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "grey15", fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand  = c(0, 0)) +
  scale_fill_manual(values = FLAG_COLS, labels = FLAG_LABS,
                    name = "Co-infection flag", drop = FALSE) +
  labs(x = NULL, y = "Proportion of BioSamples",
       title = "Co-infection flag breakdown") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

# ── Combine + save ─────────────────────────────────────────────────────────────
library(cowplot)

combined <- plot_grid(p_b, p_a, ncol = 2, align = "h", rel_widths = c(1, 1.4),
                      labels = c("A", "B"), label_size = 14)

title_row <- ggdraw() +
  draw_label("Kingdom composition of cryptic co-infections (MAL + HAL)",
             fontface = "bold", size = 15, x = 0.5, hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

final <- plot_grid(title_row, combined, ncol = 1, rel_heights = c(0.07, 1))

ggsave(paste0(OUT, ".pdf"), final, width = 14, height = 7, bg = "white",
       device = cairo_pdf)
ggsave(paste0(OUT, ".png"), final, width = 14, height = 7, bg = "white", dpi = 200)
cat("Written:", paste0(OUT, ".pdf/.png\n"))
