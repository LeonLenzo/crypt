#!/usr/bin/env Rscript
# predict_counts.R — draws the prediction-vs-annotation figure. Owns styling only.
#
# Every number comes from prep_predict_counts.py, including the axis limits, which are read
# from predict_counts_params.tsv rather than derived here. Deriving them in both places is how
# a caption ends up describing a figure that no longer exists.
#
# ASSEMBLY MODE. No panel titles, no legend, no in-panel species labels, no caption. The
# panels are printed separately and annotated in Illustrator, and predict_counts_caption.txt
# carries everything stripped from the image -- including outlier coordinates in data units.
#
# BOTH AXES OF PANELS a AND b SHARE ONE RANGE AND BREAK SET, and coord_fixed makes the panel
# square in log space. Without that the dashed 1:1 line is not at 45 degrees and the eye reads
# deviation that is not in the data, which defeats the point of a 1:1 plot.
#
#   Rscript 02_classify/report/figures/predict_counts.R

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

HERE <- "02_classify/report/figures"
d <- read.delim(file.path(HERE, "predict_counts.tsv"), stringsAsFactors = FALSE)
pr <- read.delim(file.path(HERE, "predict_counts_params.tsv"), stringsAsFactors = FALSE)
getp <- function(panel) {
  r <- pr[pr$panel == panel, ]
  list(lo = r$lo, hi = r$hi, breaks = as.numeric(strsplit(r$breaks, ";")[[1]]))
}
pa <- getp("a"); pb <- getp("b"); psz <- getp("size")

assembly_theme <- theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position  = "none",
        plot.margin      = margin(8, 12, 8, 8))

fill_scale <- scale_fill_viridis_c(option = "viridis", direction = -1, trans = "log10",
                                   limits = c(psz$lo, psz$hi), breaks = psz$breaks,
                                   oob = scales::squish)

pt <- function() geom_point(shape = 21, size = 4, stroke = 1, alpha = 0.8, colour = "black")

# One-to-one panel: shared scale on both axes, square aspect, dashed diagonal.
one_to_one <- function(p, lim, brk, labfun) {
  p + geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  linewidth = 1, colour = "grey40") +
    pt() + fill_scale +
    scale_x_log10(limits = lim, breaks = brk, labels = labfun) +
    scale_y_log10(limits = lim, breaks = brk, labels = labfun) +
    coord_fixed(ratio = 1) +
    assembly_theme
}

pA <- one_to_one(ggplot(d, aes(true_seqs, pred_seqs, fill = genome_mb)),
                 c(pa$lo, pa$hi), pa$breaks, comma) +
  labs(x = "Annotated Count", y = "Predicted Count")

pB <- one_to_one(ggplot(d, aes(true_mb, pred_mb, fill = genome_mb)),
                 c(pb$lo, pb$hi), pb$breaks, waiver()) +
  labs(x = "Annotated Mb", y = "Predicted Mb")

# Panel c is NOT 1:1 -- the two axes are different quantities -- so it keeps its own ranges
# and is left unsquared. Dotted guides at 0.90 give the eye a reference without asserting a
# threshold; 0.90 is a reading aid, not a decision rule.
pC <- ggplot(d, aes(recall, precision, fill = genome_mb)) +
  geom_hline(yintercept = 0.9, linetype = "dotted", linewidth = 0.8, colour = "grey60") +
  geom_vline(xintercept = 0.9, linetype = "dotted", linewidth = 0.8, colour = "grey60") +
  pt() + fill_scale +
  scale_x_continuous(limits = c(0.70, 1), breaks = seq(0.70, 1, 0.10)) +
  scale_y_continuous(limits = c(0.40, 1), breaks = seq(0.40, 1, 0.20)) +
  labs(x = "Recall", y = "Precision") +
  assembly_theme

p <- pA | pB | pC

ragg::agg_png(file.path(HERE, "predict_counts.png"), width = 4800, height = 1700,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "predict_counts.pdf"), width = 16, height = 5.67)
print(p); invisible(dev.off())

# Panels also written singly, since they are printed and placed separately.
for (nm in c("a", "b", "c")) {
  obj <- switch(nm, a = pA, b = pB, c = pC)
  ragg::agg_png(file.path(HERE, sprintf("predict_counts_%s.png", nm)),
                width = 1650, height = 1650, res = 300, background = "white")
  print(obj); invisible(dev.off())
  cairo_pdf(file.path(HERE, sprintf("predict_counts_%s.pdf", nm)), width = 5.5, height = 5.5)
  print(obj); invisible(dev.off())
}
cat("wrote predict_counts.{png,pdf} and predict_counts_{a,b,c}.{png,pdf}\n")
