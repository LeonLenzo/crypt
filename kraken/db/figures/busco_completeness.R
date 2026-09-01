#!/usr/bin/env Rscript
# busco_completeness.R — artistic overview of BUSCO screen results
# Output: kraken/output/figures/busco/busco_completeness.pdf + .png

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)     # geom_sina
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(scales)
})

# ── data ─────────────────────────────────────────────────────────────────────

scores <- read.delim("kraken/output/kraken_db_busco/data/busco_scores.tsv", stringsAsFactors = FALSE) |>
  filter(status %in% c("pass", "fail")) |>
  mutate(
    lineage = recode(busco_lineage,
      ascomycota_odb10    = "Ascomycota",
      basidiomycota_odb10 = "Basidiomycota",
      fungi_odb10         = "Fungi (other)",
      stramenopiles_odb10 = "Stramenopiles"
    ),
    lineage = factor(lineage,
      levels = c("Ascomycota", "Basidiomycota", "Fungi (other)", "Stramenopiles")),
    kingdom = if_else(kingdom == "oomycete", "Oomycete", "Fungal"),
    threshold = if_else(kingdom == "Oomycete", 65, 50)
  )

# per-lineage threshold lines (one value per lineage)
thresholds <- scores |>
  distinct(lineage, threshold)

n_pass  <- sum(scores$status == "pass")
n_total <- nrow(scores)

# ── palette ──────────────────────────────────────────────────────────────────

pal <- c("Fungal" = "#e8a838", "Oomycete" = "#4cb8d4")

bg      <- "#0d1117"
panel   <- "#161b22"
grid_c  <- "#30363d"
text_c  <- "#e6edf3"
subtext <- "#8b949e"

# ── plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(scores, aes(x = complete_pct, y = lineage, colour = kingdom, fill = kingdom)) +

  # --- violin background ---
  geom_violin(
    aes(group = lineage),
    width = 0.85, alpha = 0.12, linewidth = 0,
    colour = NA
  ) +

  # --- sina dots ---
  geom_sina(
    aes(group = lineage),
    size = 1.2, alpha = 0.55, stroke = 0,
    maxwidth = 0.75
  ) +

  # --- pass threshold lines ---
  geom_segment(
    data = thresholds,
    aes(x = threshold, xend = threshold,
        y = as.numeric(lineage) - 0.45,
        yend = as.numeric(lineage) + 0.45),
    colour = "#f0f6fc", linewidth = 0.55, linetype = "dashed",
    inherit.aes = FALSE
  ) +

  # --- threshold label (once per unique value) ---
  geom_text(
    data = thresholds |> distinct(threshold, .keep_all = TRUE),
    aes(x = threshold, y = lineage, label = paste0(threshold, "% threshold")),
    hjust = -0.08, vjust = -1.0,
    colour = "#f0f6fc", size = 2.8, family = "mono",
    inherit.aes = FALSE
  ) +

  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(
    limits = c(0, 102),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0, 100, 25), "%"),
    expand = expansion(mult = c(0.01, 0.05))
  ) +

  labs(
    title    = "BUSCO completeness — pathogen reference screen",
    subtitle = sprintf(
      "%d / %d assemblies pass threshold  ·  ascomycota n=1,706 markers  ·  stramenopiles n=100 markers",
      n_pass, n_total
    ),
    x = "BUSCO complete (%)",
    y = NULL
  ) +

  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = bg, colour = NA),
    panel.background  = element_rect(fill = panel, colour = NA),
    panel.grid.major.x = element_line(colour = grid_c, linewidth = 0.35),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),

    axis.text   = element_text(colour = text_c, size = 10),
    axis.title  = element_text(colour = subtext, size = 10),
    axis.ticks  = element_blank(),

    plot.title    = element_text(colour = text_c,  size = 14, face = "bold",
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(colour = subtext,  size = 8.5,
                                 margin = margin(b = 12)),
    plot.margin   = margin(16, 20, 12, 16),

    legend.position   = "bottom",
    legend.background = element_rect(fill = bg, colour = NA),
    legend.text       = element_text(colour = text_c, size = 10),
    legend.key.size   = unit(0.9, "lines"),

    strip.text = element_text(colour = text_c)
  )

# ── save ─────────────────────────────────────────────────────────────────────

out_dir <- "kraken/output/figures/busco"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "busco_completeness.pdf"),
       p, width = 9, height = 5.5, device = cairo_pdf)
ggsave(file.path(out_dir, "busco_completeness.png"),
       p, width = 9, height = 5.5, dpi = 180)

cat(sprintf("Saved to %s/busco_completeness.{pdf,png}\n", out_dir))
