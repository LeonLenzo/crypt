#!/usr/bin/env Rscript
# predict_arm_a.R — how well does ab-initio prediction recover real CDS k-mer space?
#
# Arm A of the gene-prediction benchmark: 60 annotated fungal assemblies where NCBI's own CDS
# is the ground truth, so a prediction can be scored. Written SELF-CONTAINED (axis titles,
# legend and labels kept) because it is meant to be shown to someone directly. For an
# assembly-mode version strip strip.text, axis.text.x, axis.title.y and the legend.
#
# Two panels, two different claims, and they need different displays:
#
#   A — fragmentation does NOT matter. Recall and precision are flat across four ordered
#       assembly levels. This is a NULL result, so it is drawn as boxplots over the actual
#       points. Fitting a line here would manufacture a trend out of its absence.
#
#   B — genome size DOES matter. Both metrics fall with it (r = -0.53 each, r^2 ~ 0.28), so a
#       linear fit is justified — but x is GENOME SIZE from the downloaded FASTA, not the
#       over-prediction ratio. Precision is |T n P| / |P|, so plotting it against anything
#       derived from |P| is close to tautological. Genome size is independent of both sketches
#       and, in fungi, is mostly a proxy for transposon content.
#
# The mechanism the outliers show: ab-initio prediction calls genes inside transposons.
# Blumeria graminis is ~64% repetitive and was over-predicted 3.91x.
#
# Note r^2 ~ 0.28 — genome size explains under a third of the variance, so the fit indicates
# direction and should not be read as prediction. The scatter is deliberately drawn on top.
#
#   Rscript 02_classify/report/figures/predict_arm_a.R

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel)
})

HERE <- "02_classify/report/figures"
# base read.delim, not readr: vroom.so fails to load on this machine
d <- read.delim(file.path(HERE, "predict_arm_a.tsv"), stringsAsFactors = FALSE)

# "Complete Genome" collides with "Chromosome" on a shared axis at this width.
LEVELS <- c("Contig", "Scaffold", "Chromosome", "Complete")
d$assembly_level[d$assembly_level == "Complete Genome"] <- "Complete"
d$assembly_level <- factor(d$assembly_level, levels = LEVELS)

# Colourblind-safe default pair. Colour means METRIC and is consistent across both panels.
COL <- c(recall = "#2a78d6", precision = "#eb6834")
# Legend keys stay SHORT. The full definitions ran off the canvas at base_size 16; they
# belong in the caption, not in a key.
LAB <- c(recall = "Recall", precision = "Precision")

long <- rbind(
  data.frame(d[c("genome_mb", "assembly_level", "organism_name", "focal")],
             metric = "recall", value = d$recall),
  data.frame(d[c("genome_mb", "assembly_level", "organism_name", "focal")],
             metric = "precision", value = d$precision))
long$metric <- factor(long$metric, levels = c("recall", "precision"))

base_theme <- theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position  = "top",
        legend.title     = element_blank(),
        plot.subtitle    = element_text(size = 14),
        plot.margin      = margin(6, 14, 6, 6))

# ── Panel A: flat across fragmentation ──────────────────────────────────────────────────────
pA <- ggplot(long, aes(assembly_level, value, fill = metric)) +
  geom_boxplot(alpha = 0.3, fill = "lightgray", colour = "black",
               linewidth = 1, width = 0.6, outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  geom_point(aes(fill = metric), shape = 21, size = 4, stroke = 1, alpha = 0.7,
             colour = "black",
             position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75)) +
  scale_fill_manual(values = COL, labels = LAB) +
  scale_y_continuous(breaks = seq(0.4, 1, 0.1)) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(x = NULL, y = "Containment against true CDS",
       subtitle = "Fragmentation does not degrade recovery") +
  base_theme

# ── Panel B: both metrics fall with genome size ─────────────────────────────────────────────
# Label the cases that drive the trend, plus anything with precision under 0.8.
lab <- long[long$metric == "precision" & (long$value < 0.82 | long$genome_mb > 100), ]
italicise <- function(x) sprintf('italic("%s")', x)   # plotmath; ggtext dies under theme_minimal

pB <- ggplot(long, aes(genome_mb, value, fill = metric)) +
  # show.legend = FALSE: geom_smooth inherits fill for its ribbon, which puts a line glyph in
  # the fill key. pA's key has points only, the two guide sets then differ, and
  # plot_layout(guides = "collect") gives up and draws the legend twice.
  geom_smooth(aes(colour = metric), method = "lm", formula = y ~ x,
              se = TRUE, linewidth = 1.2, alpha = 0.15, show.legend = FALSE) +
  geom_point(shape = 21, size = 4, stroke = 1, alpha = 0.7, colour = "black") +
  geom_text_repel(data = lab, aes(label = italicise(organism_name)),
                  parse = TRUE, size = 3.6, colour = "grey25",
                  min.segment.length = 0, segment.colour = "grey60",
                  box.padding = 0.8, point.padding = 0.4, force = 3,
                  max.overlaps = 20, seed = 1) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160)) +
  scale_fill_manual(values = COL, labels = LAB) +
  scale_colour_manual(values = COL, labels = LAB) +
  # pA guides fill only; without this pB also guides colour, the guide sets differ and
  # plot_layout(guides = "collect") draws the legend TWICE instead of merging it.
  guides(colour = "none") +
  scale_y_continuous(breaks = seq(0.4, 1, 0.1)) +
  # coord_cartesian, NOT scale limits: the lm confidence ribbon for recall rises above 1.0 at
  # small genome sizes, and a scale limit DELETES those rows (warning: "Removed 32 rows") and
  # with them part of the ribbon. coord_cartesian zooms the view and fits on all 60 points.
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(x = "Genome size (Mb, log scale)", y = "Containment against true CDS",
       subtitle = expression("Larger genomes are harder, "*r^2*" = 0.28")) +
  base_theme

p <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "60 annotated fungal assemblies; NCBI CDS is ground truth. FracMinHash k=31, scaled=200.",
      "Recall = fraction of true CDS k-mer space recovered; precision = fraction of predicted",
      "k-mer space that is real CDS. Labelled: precision < 0.82 or genome > 100 Mb.", sep = "\n"),
    theme = theme(plot.caption = element_text(size = 11, hjust = 0, colour = "grey30"))) &
  theme(plot.tag = element_text(size = 18, face = "bold"), legend.position = "top")

ragg::agg_png(file.path(HERE, "predict_arm_a.png"), width = 3600, height = 1600,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "predict_arm_a.pdf"), width = 12, height = 5.33)
print(p); invisible(dev.off())
cat("wrote predict_arm_a.png and .pdf\n")
