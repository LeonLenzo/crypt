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
pa <- getp("a"); pb <- getp("b"); pc <- getp("c"); psz <- getp("size")
A_UNIT <- getp("a_unit")$lo          # panel a is drawn in units of 10,000

# theme_classic, not theme_minimal: real black axis LINES with ticks and no grid. The house
# theme is minimal, but Leon asked for drawn axes here and a grid behind a 1:1 plot competes
# with the diagonal, which is the line that matters.
assembly_theme <- theme_classic(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 1),
        axis.ticks       = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(4, "pt"),
        # Tick numerals a little above base_size: they are the only text surviving into the
        # assembled figure, so they carry the reading on their own.
        axis.text        = element_text(colour = "black", size = 18),
        legend.position  = "none",
        plot.margin      = margin(8, 12, 8, 8))

fill_scale <- scale_fill_viridis_c(option = "viridis", direction = -1, trans = "log10",
                                   limits = c(psz$lo, psz$hi), breaks = psz$breaks,
                                   oob = scales::squish)

pt <- function() geom_point(shape = 21, size = 4, stroke = 1, alpha = 0.8, colour = "black")

# Pull the guide grob out of a built plot. ggplot2 4.x renames the legend grob, so match on
# either spelling rather than assuming one.
cowplot_get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  i <- which(vapply(g$grobs, function(x) !is.null(x$name) &&
                      grepl("guide-box|legend", x$name), logical(1)))
  if (!length(i)) stop("no legend grob found")
  g$grobs[[i[1]]]
}

# One-to-one panel: shared ZERO-BASED linear scale on both axes, square, dashed diagonal.
# expand = c(0, 0) so the axes meet at the origin rather than floating off it.
one_to_one <- function(p, lim, brk, labfun) {
  p + geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  linewidth = 1, colour = "grey40") +
    pt() + fill_scale +
    scale_x_continuous(limits = lim, breaks = brk, labels = labfun,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = lim, breaks = brk, labels = labfun,
                       expand = expansion(mult = c(0, 0.02))) +
    coord_fixed(ratio = 1) +
    assembly_theme
}

# Numerals divided by 10,000 and the factor moved into the axis title, so the ticks read
# 0, 0.5, 1.0 ... instead of 0, 5,000, 10,000 ... which collided at this panel width.
in_units <- function(x) ifelse(is.na(x), NA, format(x / A_UNIT, trim = TRUE))
pA <- one_to_one(ggplot(d, aes(true_seqs, pred_seqs, fill = genome_mb)),
                 c(pa$lo, pa$hi), pa$breaks, in_units) +
  # Plain text, not plotmath: expression(...%*%...) renders as a CENTRE DOT here, not a
  # multiplication sign. ragg handles the unicode glyph directly.
  labs(x = "Annotated Count (\u00d7 10,000)", y = "Predicted Count (\u00d7 10,000)")

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
  # ONE range on both axes so recall and precision are read on the same scale, and squared
  # for the same reason. It reaches to 0.4 because the lowest precision is 0.487, well below
  # any recall.
  scale_x_continuous(limits = c(pc$lo, pc$hi), breaks = pc$breaks,
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(pc$lo, pc$hi), breaks = pc$breaks,
                     expand = expansion(mult = 0.02)) +
  coord_fixed(ratio = 1) +
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
# The legend is stripped from every panel, so it is written separately for placement. Vertical
# as asked; the figure is three panels wide and a vertical bar sits beside them without
# stealing width.
legend_only <- ggplot(d, aes(true_mb, pred_mb, fill = genome_mb)) +
  geom_point(shape = 21, size = 4, stroke = 1, colour = "black") +
  scale_fill_viridis_c(name = "Genome\nsize (Mb)", option = "viridis", direction = -1,
                       trans = "log10", limits = c(psz$lo, psz$hi), breaks = psz$breaks,
                       labels = psz$breaks, oob = scales::squish,
                       guide = guide_colourbar(barwidth = 0.9, barheight = 11,
                                               frame.colour = "black", frame.linewidth = 0.6,
                                               ticks.colour = "black", title.position = "top")) +
  theme_classic(base_size = 16) +
  theme(legend.position = "right", legend.background = element_blank(),
        legend.text = element_text(colour = "black"))
leg <- cowplot_get_legend(legend_only)
ragg::agg_png(file.path(HERE, "predict_counts_legend.png"), width = 700, height = 1400,
              res = 300, background = "white")
grid::grid.newpage(); grid::grid.draw(leg); invisible(dev.off())
cairo_pdf(file.path(HERE, "predict_counts_legend.pdf"), width = 2.33, height = 4.67)
grid::grid.newpage(); grid::grid.draw(leg); invisible(dev.off())

cat("wrote predict_counts.{png,pdf}, predict_counts_{a,b,c}.{png,pdf}",
    "and predict_counts_legend.{png,pdf}\n")
