#!/usr/bin/env Rscript
# Kraken2 pilot: three-way comparison STAT / kallisto / Kraken2
# Panel A: STAT euk_pct vs Kraken2  — shows the STAT blind spot
# Panel B: kallisto vs Kraken2      — validates Kraken2 against kallisto
# Run from crypt/: Rscript kraken/pilot/kraken/figure_kraken.R

library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

dat <- read.delim("kraken/pilot/kraken/kraken_results.tsv", stringsAsFactors = FALSE) |>
  mutate(
    # total classified by Kraken2 (any of the 3 target organisms)
    # PGT reads cross-classify to PST (same-genus k-mer sharing), so
    # pct_classified is more appropriate than the organism-specific column
    kraken_pct = pct_classified,
    organism = case_when(
      grepl("striiformis", organism) ~ "P. striiformis (PST)",
      grepl("graminis",    organism) ~ "P. graminis (PGT)",
      grepl("oryzae",      organism) ~ "P. oryzae (POR)"
    ),
    tier = factor(tier, levels = c("zero", "low", "high"),
                  labels = c("Zero (STAT = 0%)", "Low (STAT 0-1%)", "High (STAT >= 1%)")),
    # label the two extreme PST blind-spot runs in panel A
    label_a = ifelse(organism == "P. striiformis (PST)" &
                       tier == "Zero (STAT = 0%)" & kraken_pct > 20, Run, NA),
    # label PGT high-tier in panel B (ref quality outliers)
    label_b = ifelse(organism == "P. graminis (PGT)" &
                       tier == "High (STAT >= 1%)" & kallisto_pct > 10, Run, NA)
  )

org_colours <- c(
  "P. striiformis (PST)" = "#e67e22",
  "P. graminis (PGT)"    = "#8e44ad",
  "P. oryzae (POR)"      = "#27ae60"
)

tier_shapes <- c(
  "Zero (STAT = 0%)"  = 21,
  "Low (STAT 0-1%)"   = 24,
  "High (STAT >= 1%)" = 22
)

log_scale <- list(
  scale_x_continuous(trans   = scales::log1p_trans(),
                     limits  = c(0, 110),
                     breaks  = c(1, 10, 100),
                     labels  = c("0", "1", "2"),
                     expand  = expansion(mult = c(0.01, 0.03))),
  scale_y_continuous(trans   = scales::log1p_trans(),
                     limits  = c(0, 110),
                     breaks  = c(1, 10, 100),
                     labels  = c("0", "1", "2"),
                     expand  = expansion(mult = c(0.01, 0.03))),
  coord_cartesian(xlim = c(0, 105), ylim = c(0, 105))
)

base_theme <- list(
  scale_fill_manual(values = org_colours, name = "Organism"),
  scale_shape_manual(values = tier_shapes, name = "STAT tier"),
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(fill  = "grey50", size = 4))),
  theme_classic(base_size = 11),
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.9, "lines"),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    plot.tag         = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3)
  )
)

# ── Panel A: STAT vs kallisto ─────────────────────────────────────────────────

pA <- ggplot(dat, aes(x = euk_pct, y = kallisto_pct,
                      fill = organism, shape = tier)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.5) +
  geom_point(size = 3.5, colour = "white", stroke = 0.4, alpha = 0.9) +
  geom_label_repel(aes(label = label_a), fill = "white", size = 2.4,
                   colour = "#e67e22", label.padding = 0.15,
                   box.padding = 0.4, max.overlaps = 20,
                   na.rm = TRUE, show.legend = FALSE) +
  log_scale +
  base_theme +
  labs(
    tag   = "A",
    x     = expression("STAT eukaryotic pathogen, log"[10] * "(%)"),
    y     = expression("Kallisto pseudoalignment rate, log"[10] * "(%)"),
    title = "STAT vs kallisto",
    subtitle = "PST runs sit at STAT = 0% regardless of true abundance"
  )

# ── Panel B: kallisto vs Kraken2 ─────────────────────────────────────────────

pB <- ggplot(dat, aes(x = kallisto_pct, y = kraken_pct,
                      fill = organism, shape = tier)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.5) +
  geom_point(size = 3.5, colour = "white", stroke = 0.4, alpha = 0.9) +
  geom_label_repel(aes(label = label_b), fill = "white", size = 2.4,
                   colour = "#8e44ad", label.padding = 0.15,
                   box.padding = 0.4, max.overlaps = 20,
                   na.rm = TRUE, show.legend = FALSE) +
  log_scale +
  base_theme +
  labs(
    tag   = "B",
    x     = expression("Kallisto pseudoalignment rate, log"[10] * "(%)"),
    y     = expression("Kraken2 target organism, log"[10] * "(%)"),
    title = "Kallisto vs Kraken2",
    subtitle = "PGT outliers (purple, annotated): same-genus k-mer sharing with PST"
  )

# ── Combine ───────────────────────────────────────────────────────────────────

combined <- pA + pB +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Three-method comparison: STAT / kallisto / Kraken2  |  45 MAL pilot runs",
    subtitle = "500k reads/run streamed from ENA  ·  3 rust/blast organisms  ·  3 STAT signal tiers",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

ggsave("kraken/pilot/kraken/figure_kraken.pdf", combined, width = 12, height = 5)
ggsave("kraken/pilot/kraken/figure_kraken.png", combined, width = 12, height = 5, dpi = 200)
cat("Written: kraken/pilot/kraken/figure_kraken.pdf/.png\n")
