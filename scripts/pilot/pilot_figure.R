#!/usr/bin/env Rscript
# Kallisto pilot: STAT euk_pct vs kallisto pseudoalignment rate
# 45 runs × 3 organisms × 3 tiers — diagnostic record of STAT reliability
# Run from crypt/: Rscript scripts/pilot/pilot_figure.R

library(ggplot2)
library(dplyr)
library(ggrepel)

dat <- read.delim("scripts/pilot/results.tsv", stringsAsFactors = FALSE) |>
  mutate(
    organism = sub("Puccinia striiformis f. sp. tritici", "P. striiformis (PST)", organism),
    organism = sub("Puccinia graminis f. sp. tritici",   "P. graminis (PGT)",    organism),
    organism = sub("Pyricularia oryzae",                 "P. oryzae (POR)",      organism),
    tier = factor(tier, levels = c("zero", "low", "high"),
                  labels = c("Zero (STAT = 0%)", "Low (STAT 0-1%)", "High (STAT >= 1%)")),
    acc_type = ifelse(grepl("^ERR|^DRR", Run), "ERR/DRR", "SRR"),
    # flag the notable PST outliers for annotation
    label = ifelse(organism == "P. striiformis (PST)" &
                     tier != "High (STAT >= 1%)" & kallisto_pct > 20, Run, NA)
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

p <- ggplot(dat, aes(x = euk_pct, y = kallisto_pct,
                     fill = organism, shape = tier)) +
  # y = x reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.5) +
  # points
  geom_point(size = 3.5, colour = "white", stroke = 0.4, alpha = 0.9) +
  # annotate PST outliers
  geom_label_repel(aes(label = label), fill = "white", size = 2.5,
                   colour = "#e67e22", label.padding = 0.15,
                   box.padding = 0.4, max.overlaps = 20,
                   na.rm = TRUE, show.legend = FALSE) +
  scale_fill_manual(values = org_colours, name = "Organism") +
  scale_shape_manual(values = tier_shapes, name = "STAT tier") +
  scale_x_continuous(trans = scales::log1p_trans(),
                     limits = c(0, 110),
                     breaks = c(1, 10, 100),
                     labels = c("0", "1", "2"),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(trans = scales::log1p_trans(),
                     limits = c(0, 110),
                     breaks = c(1, 10, 100),
                     labels = c("0", "1", "2"),
                     expand = expansion(mult = c(0.01, 0.03))) +
  coord_cartesian(xlim = c(0, 105), ylim = c(0, 105)) +
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(fill  = "grey50", size = 4))) +
  labs(
    x     = expression("STAT eukaryotic pathogen reads, log"[10] * "(%)"),
    y     = expression("Kallisto pseudoalignment rate, log"[10] * "(%)"),
    title = "STAT vs kallisto: pathogen detection in MAL runs",
    subtitle = paste0(
      "45 runs x 3 organisms x 3 STAT signal tiers  |  100k reads/run\n",
      "Dashed line = perfect agreement. PST outlier runs annotated (all ERR accessions)."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.9, "lines"),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3)
  )

ggsave("scripts/pilot/pilot_figure.pdf", p, width = 7, height = 5.5)
ggsave("scripts/pilot/pilot_figure.png", p, width = 7, height = 5.5, dpi = 200)
cat("Written: scripts/pilot/pilot_figure.pdf/.png\n")
