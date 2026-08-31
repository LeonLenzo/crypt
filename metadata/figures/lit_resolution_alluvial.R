#!/usr/bin/env Rscript
# metadata/figures/lit_resolution_alluvial.R — literature resolution funnel (ggalluvial)
#
# Reads metadata/output/figures/sankey/lit_resolution_data.tsv (built by
# prep_lit_resolution.py) and plots the meta_search.py -> meta_text.py
# resolution funnel as an alluvial diagram: three axes (how the DOI was
# found, how full text was retrieved, final outcome), one ribbon per
# BioProject, coloured by final outcome.
#
# Full manual control over layout (unlike the retired Plotly version):
# node order = factor levels below, node/ribbon positions are real ggplot2
# coordinates. Output is PNG + SVG (SVG is hand-editable in Illustrator/
# Inkscape if you want to nudge things further — no interactive HTML, but
# no fighting an auto-layout engine either).
#
# Run from crypt/: Rscript metadata/figures/lit_resolution_alluvial.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
})

IN_TSV  <- "metadata/output/figures/sankey/lit_resolution_data.tsv"
OUT_PNG <- "metadata/output/figures/sankey/lit_resolution_alluvial.png"
OUT_SVG <- "metadata/output/figures/sankey/lit_resolution_alluvial.svg"

d <- read.delim(IN_TSV, stringsAsFactors = FALSE)

# Node order per axis — edit these to reorder nodes top-to-bottom
search_levels  <- c("NCBI (PMID)", "NCBI (DOI only)", "Web search", "No DOI found")
text_levels    <- c("PMC OA", "Unpaywall PDF", "Manual PDF", "No OA copy found", "No DOI")
outcome_levels <- c("Full text retrieved", "No full text available")

d$search  <- factor(d$search,  levels = search_levels)
d$text    <- factor(d$text,    levels = text_levels)
d$outcome <- factor(d$outcome, levels = outcome_levels)

n_total <- nrow(d)
n_full  <- sum(d$outcome == "Full text retrieved")
n_none  <- sum(d$outcome == "No full text available")

INCLUDED <- "#2980b9"
EXCLUDED <- "#95a5a6"

p <- ggplot(d, aes(axis1 = search, axis2 = text, axis3 = outcome)) +
  geom_alluvium(aes(fill = outcome), width = 1/8, alpha = 0.7, knot.pos = 0.4) +
  geom_stratum(width = 1/8, fill = "grey93", color = "white") +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(paste0(stratum, "\n", count))),
    size = 3.3, lineheight = 0.85
  ) +
  scale_x_discrete(limits = c("Search", "Text", "Outcome"), expand = c(.1, .1)) +
  scale_fill_manual(values = c(
    "Full text retrieved"    = INCLUDED,
    "No full text available" = EXCLUDED
  )) +
  labs(
    title = sprintf(
      "Literature resolution — meta_search.py → meta_text.py — %s BioProjects",
      format(n_total, big.mark = ",")
    ),
    subtitle = sprintf(
      "Full text retrieved: %d (%.1f%%)   |   No full text available: %d (%.1f%%)",
      n_full, 100 * n_full / n_total, n_none, 100 * n_none / n_total
    ),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
  )

ggsave(OUT_PNG, p, width = 12, height = 7, dpi = 150, bg = "white")
ggsave(OUT_SVG, p, width = 12, height = 7, bg = "white")
cat("Written", OUT_PNG, "\n")
cat("Written", OUT_SVG, "\n")
