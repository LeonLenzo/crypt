#!/usr/bin/env Rscript
# predict_counts.R — does ab-initio prediction reproduce NCBI's annotation?
#
# Three panels asking three increasingly strict questions of the same 60 assemblies:
#
#   a  Does it call the same NUMBER of CDS?
#   b  Does it call the same AMOUNT of sequence?
#   c  Is it the SAME sequence?
#
# a and b can both agree while c fails completely — identical totals say nothing about
# identity, which is the whole reason the benchmark measures k-mer containment rather than
# stopping at counts. Panels a and b carry a dashed 1:1 line, so deviation is read directly
# off the diagonal rather than inferred from a fitted slope.
#
# Counts are fair to compare: NCBI's cds_from_genomic.fna holds one entry per CDS INCLUDING
# isoforms, while vipsania writes one per predicted gene. Checked rather than assumed — the
# ratio of NCBI CDS entries to its own protein-coding gene count has median 1.000 across the
# 60, and only 4 exceed 1.10. Fungal isoform inflation is negligible at this scale.
#
# Colour is genome size throughout, because that is what the earlier analysis showed drives
# the deviation (r = -0.53 against both recall and precision): ab-initio prediction calls
# genes inside transposons, and large fungal genomes are large mostly because of transposons.
#
#   Rscript 02_classify/report/figures/predict_counts.R

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

HERE <- "02_classify/report/figures"
d <- read.delim(file.path(HERE, "predict_counts.tsv"), stringsAsFactors = FALSE)

d$ratio_seqs  <- d$pred_seqs / d$true_seqs
d$ratio_bases <- d$pred_mb   / d$true_mb
italicise <- function(x) sprintf('italic("%s")', x)   # plotmath; ggtext dies under theme_minimal

base_theme <- theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        plot.subtitle    = element_text(size = 13),
        legend.position  = "top",
        plot.margin      = margin(6, 14, 6, 6))

# Continuous fill for genome size. Marks stay house-style: shape 21, black border, size 4,
# stroke 1, alpha 0.7.
fill_scale <- scale_fill_viridis_c(
  name = "Genome size (Mb)", option = "viridis", direction = -1,
  trans = "log10", breaks = c(10, 20, 40, 80, 160),
  guide = guide_colourbar(barwidth = 9, barheight = 0.7, title.vjust = 0.9))

pt <- function() geom_point(shape = 21, size = 4, stroke = 1, alpha = 0.8, colour = "black")
diag_line <- function() geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                                    linewidth = 1, colour = "grey40")

# ── a: number of CDS called ─────────────────────────────────────────────────────────────────
lab_a <- d[d$ratio_seqs > 1.6 | d$ratio_seqs < 0.6, ]
pA <- ggplot(d, aes(true_seqs, pred_seqs, fill = genome_mb)) +
  diag_line() + pt() +
  geom_text_repel(data = lab_a, aes(label = italicise(organism_name)), parse = TRUE,
                  size = 3.4, colour = "grey25", min.segment.length = 0,
                  segment.colour = "grey60", box.padding = 0.9, point.padding = 0.4, force = 4,
                  max.overlaps = 20, seed = 1) +
  scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
  fill_scale +
  labs(x = "CDS called by NCBI annotation", y = "CDS called by prediction",
       subtitle = "a  Number of CDS: median ratio 1.01") +
  base_theme

# ── b: amount of sequence ───────────────────────────────────────────────────────────────────
lab_b <- d[d$ratio_bases > 1.45 | d$ratio_bases < 0.75, ]
pB <- ggplot(d, aes(true_mb, pred_mb, fill = genome_mb)) +
  diag_line() + pt() +
  geom_text_repel(data = lab_b, aes(label = italicise(organism_name)), parse = TRUE,
                  size = 3.4, colour = "grey25", min.segment.length = 0,
                  segment.colour = "grey60", box.padding = 0.9, point.padding = 0.4, force = 4,
                  max.overlaps = 20, seed = 1) +
  scale_x_log10() + scale_y_log10() +
  fill_scale +
  labs(x = "CDS Mb from NCBI annotation", y = "CDS Mb from prediction",
       subtitle = "b  Amount of sequence: median ratio 1.03") +
  base_theme

# ── c: is it the SAME sequence? ─────────────────────────────────────────────────────────────
# The strict question. Equal totals in a and b are consistent with completely different
# sequence; only containment distinguishes them. Top-right is agreement.
d$recall <- as.numeric(d$recall); d$precision <- as.numeric(d$precision)
lab_c <- d[d$precision < 0.80 | d$recall < 0.87, ]
pC <- ggplot(d, aes(recall, precision, fill = genome_mb)) +
  geom_hline(yintercept = 0.9, linetype = "dotted", linewidth = 0.8, colour = "grey60") +
  geom_vline(xintercept = 0.9, linetype = "dotted", linewidth = 0.8, colour = "grey60") +
  pt() +
  geom_text_repel(data = lab_c, aes(label = italicise(organism_name)), parse = TRUE,
                  size = 3.4, colour = "grey25", min.segment.length = 0,
                  segment.colour = "grey60", box.padding = 0.9, point.padding = 0.4, force = 4,
                  max.overlaps = 20, seed = 1) +
  scale_x_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.1)) +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  fill_scale +
  # Short axis titles: the full definitions are longer than the panel is tall and get clipped
  # when rotated. They live in the caption instead.
  labs(x = "Recall", y = "Precision",
       subtitle = "c  Same sequence? median 0.96 / 0.94") +
  base_theme

p <- (pA | pB | pC) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = paste(
      "60 annotated fungal assemblies; NCBI CDS is ground truth. Dashed line in a and b is 1:1.",
      "c: recall = fraction of true CDS k-mer space recovered; precision = fraction of predicted k-mer space that is real CDS.",
      "Containment from FracMinHash sketches, k=31, scaled=200.",
      "Equal totals in a and b would be consistent with entirely different sequence; only c tests identity.",
      sep = "\n"),
    theme = theme(plot.caption = element_text(size = 11, hjust = 0, colour = "grey30"))) &
  # bottom, not top: collected at "top" the colourbar centres over panel b and crowds its
  # plot area, while the outer panels start higher and look unbalanced.
  theme(legend.position = "bottom", legend.box.margin = margin(t = -2))

ragg::agg_png(file.path(HERE, "predict_counts.png"), width = 5100, height = 1750,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "predict_counts.pdf"), width = 17, height = 5.83)
print(p); invisible(dev.off())
cat("wrote predict_counts.png and .pdf\n")
