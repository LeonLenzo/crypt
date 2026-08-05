#!/usr/bin/env Rscript
# figure/scatter/scatter.R — MAL and HAL scatter panels side by side
# Run from crypt/: Rscript figure/scatter/scatter.R
# Requires scatter_data.tsv from prep_scatter.py

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

DAT <- "figure/scatter/scatter_data.tsv"
OUT <- "figure/scatter/scatter"

LOG_FLOOR <- 0.05

MAL_COLS <- c("MAL_single" = "#7fb3d3", "MAL_coinf" = "#1a5276")
HAL_COLS <- c("HAL_single" = "#82c995", "HAL_coinf" = "#166420")

MAL_LABS <- c("MAL_single" = "Single pathogen", "MAL_coinf" = "Co-infected")
HAL_LABS <- c("HAL_single" = "Single pathogen", "HAL_coinf" = "Co-infected")

# ── Load data ─────────────────────────────────────────────────────────────────
d <- read.delim(DAT, stringsAsFactors = FALSE) |>
  mutate(
    host_pct = pmax(host_pct, LOG_FLOOR),
    euk_pct  = pmax(euk_pct,  LOG_FLOOR),
    group    = case_when(
      layer == "gate" ~ paste0(mode, "_", ifelse(coinf == "True", "coinf", "single")),
      TRUE            ~ "grey"
    )
  )

grey_mal <- filter(d, group == "grey", mode == "MAL")
grey_hal <- filter(d, group == "grey", mode == "HAL")
mal_pts  <- filter(d, group %in% names(MAL_COLS)) |>
              mutate(group = factor(group, levels = names(MAL_COLS))) |>
              arrange(group)
hal_pts  <- filter(d, group %in% names(HAL_COLS)) |>
              mutate(group = factor(group, levels = names(HAL_COLS))) |>
              arrange(group)

cat("MAL:", nrow(mal_pts), " grey:", nrow(grey_mal),
    " HAL:", nrow(hal_pts), " grey:", nrow(grey_hal), "\n")
cat("MAL breakdown:\n"); print(as.data.frame(count(mal_pts, group)))
cat("HAL breakdown:\n"); print(as.data.frame(count(hal_pts, group)))

# ── Shared axis + theme ───────────────────────────────────────────────────────
shared_x <- scale_x_log10(
  breaks = c(0.1, 1, 10, 100),
  labels = c("0.1%", "1%", "10%", "100%"),
  limits = c(LOG_FLOOR, 100)
)
shared_y <- scale_y_log10(
  breaks = c(0.1, 1, 10, 100),
  labels = c("0.1%", "1%", "10%", "100%"),
  limits = c(LOG_FLOOR, 100)
)
base_theme <- theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

grey_layer_mal <- geom_point(
  data = grey_mal, aes(x = host_pct, y = euk_pct),
  size = 1, alpha = 0.12, colour = "#8c9a9e", stroke = 0
)
grey_layer_hal <- geom_point(
  data = grey_hal, aes(x = host_pct, y = euk_pct),
  size = 1, alpha = 0.12, colour = "#8c9a9e", stroke = 0
)

# ── MAL panel: gate line at x = 1 (Viridiplantae >= 1%) ──────────────────────
p_mal <- ggplot() +
  grey_layer_mal +
  geom_point(data = mal_pts,
             aes(x = host_pct, y = euk_pct, colour = group),
             size = 1.0, alpha = 0.6, stroke = 0) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  shared_x + shared_y +
  scale_colour_manual(values = MAL_COLS, labels = MAL_LABS, name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(
    title    = "MAL — microbe-as-library",
    subtitle = sprintf("Gate: Host ≥ 1%%  ·  n = %s BioSamples",
                       formatC(nrow(mal_pts), format = "d", big.mark = ",")),
    x = "Host k-mers (%)",
    y = "Pathogen k-mers (%)"
  ) +
  base_theme

# ── HAL panel: gate line at y = 1 (euk pathogen >= 1%) ───────────────────────
p_hal <- ggplot() +
  grey_layer_hal +
  geom_point(data = hal_pts,
             aes(x = host_pct, y = euk_pct, colour = group),
             size = 1.0, alpha = 0.6, stroke = 0) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  shared_x + shared_y +
  scale_colour_manual(values = HAL_COLS, labels = HAL_LABS, name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(
    title    = "HAL — host-as-library",
    subtitle = sprintf("Gate: Pathogen ≥ 1%%  ·  n = %s BioSamples",
                       formatC(nrow(hal_pts), format = "d", big.mark = ",")),
    x = "Host k-mers (%)",
    y = NULL
  ) +
  base_theme

# ── Combine ───────────────────────────────────────────────────────────────────
p <- p_mal + p_hal +
  plot_annotation(
    title    = "STAT defined Host vs Pathogen k-mers",
    subtitle = sprintf(
      "Grey: mode-specific failed-gate runs  ·  Coloured: %s gate-pass BioSamples",
      formatC(nrow(mal_pts) + nrow(hal_pts), format = "d", big.mark = ",")
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(colour = "grey50", size = 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave(paste0(OUT, ".pdf"), p, width = 14, height = 7, bg = "white",
       device = cairo_pdf)
ggsave(paste0(OUT, ".png"), p, width = 14, height = 7, bg = "white", dpi = 200)
cat("Written:", paste0(OUT, ".pdf/.png\n"))
