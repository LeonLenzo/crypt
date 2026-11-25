#!/usr/bin/env Rscript
# plot_runs.R — four columns of samples.png, opened up on the axes the call is made on.
#
#   python release/03_classify/figures/runs/prep_runs.py   # on Setonix; decides
#   Rscript release/03_classify/figures/runs/plot_runs.R   # draws
#
# Reads  runs_points.tsv, runs_floors.tsv, ../../data/inputs/panel_runs.tsv
# Writes runs_<rung>_<label>.png, one per panel, and runs_legend.png
#
# Same split as the samples panel: prep decides, this draws. Every colour arrives as a hex
# computed from 02_reference's MATRIX, so a category is one colour across the chapter.
#
# ## The axes, which match the research-tree figure this replaces
#
# Redundancy available against coverage observed, both normalised by the taxon's own db2
# minimizer total, which is what makes the null ONE curve rather than one per taxon:
#
#   solid  y = 1-exp(-x)            what random sampling at that depth would give
#   dashed y = 0.05*(1-exp(-x))     the depth-normalised abundance floor
#   dotted y = 0.02                 the coverage floor, anchored on the off-host null
#
# Drawing the saturation floor as a CURVE and not a horizontal line is the point of the
# 2026-11-24 measurement: it is an abundance floor that bends down at low depth, and the two
# floors cross at r = 0.51. Above that the dashed curve binds, below it the dotted line does.
#
# ## Assembly mode, and SEPARATE FILES
#
# leon-figure-style: black axes and ticks, no grid, no panel titles, no taxon labels. The
# AXIS titles stay, as they are on the panels in research/02_classify these replace. These are
# placed down the side of samples.png and captioned there, so a title on each would be read
# four times and the species names would be covered by the main panel. One file per panel so
# they can be placed independently. Shared limits across all four, because the comparison
# between rungs is the point and four self-scaled panels would hide it.

suppressPackageStartupMessages({library(ggplot2); library(grid); library(scales)})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(HERE) || HERE == "") HERE <- "release/03_classify/figures/runs"
rd <- function(f, ...) read.delim(file.path(HERE, f), sep = "\t", quote = "",
                                  check.names = FALSE, colClasses = "character", ...)

pts <- rd("runs_points.tsv")
fl  <- rd("runs_floors.tsv")
pan <- read.delim(file.path(HERE, "../../data/inputs/panel_runs.tsv"), sep = "\t",
                  comment.char = "#", colClasses = "character")

pts$redundancy <- as.numeric(pts$redundancy)
pts$coverage   <- as.numeric(pts$coverage)
pts$rung       <- as.integer(pts$rung)
pts$called     <- pts$call == "detected"
gv <- function(k) as.numeric(fl$value[fl$item == k])
CF <- gv("coverage_floor"); SF <- gv("saturation_floor")

null <- fl[grepl("^null:", fl$item), ]
null <- data.frame(x = as.numeric(sub("^null:", "", null$item)),
                   y = as.numeric(null$value))

pan$rung <- as.integer(pan$rung)
pan <- pan[order(pan$rung), ]
fills <- sort(unique(pts$fill)); edges <- sort(unique(pts$edge))
# The y limit runs to 1, not to the largest observed coverage. Stopping at the data max
# clipped the NULL CURVE, which rises to 100% by definition: the curve is the reference the
# points are read against, so cutting it off is cutting off the thing being compared to.
XL <- range(pts$redundancy)
YL <- c(min(pts$coverage), 1)

# `which_rung`, NOT `rung`: subset() evaluates its expression inside the data frame, so a
# parameter sharing a column's name resolves to the COLUMN and the filter is always TRUE.
# Every panel came out holding every run's points before this was caught.
panel_plot <- function(which_rung) {
  d <- pts[pts$rung == which_rung, ]
  ggplot(d, aes(redundancy, coverage)) +
    # THE NULL CURVE IS NOT DRAWN, Leon's call. y = 1-exp(-x) is what sampling at that
    # depth would give if every read landed on a new k-mer, and no point reaches it -- real
    # transcriptomes are not uniform over a genome's k-mer space -- so as a line on this
    # chart it marked a ceiling nothing approaches and competed with the two floors that
    # DO decide calls. It is still the thing the abundance floor is a fraction OF, which
    # is why that floor is a curve and not a line, and the caption says so.
    geom_line(data = null, aes(x, y * SF), inherit.aes = FALSE,
              colour = "#4a4a4a", linewidth = 0.5, linetype = "22") +
    geom_hline(yintercept = CF, colour = "#4a4a4a", linewidth = 0.4, linetype = "12") +
    geom_point(data = d[!d$called, ], shape = 21, size = 0.8, stroke = 0.2,
               fill = "#ffffff", colour = "#c4c4c4") +
    geom_point(data = d[d$called, ], aes(fill = fill, colour = edge),
               shape = 21, size = 2.6, stroke = 0.5) +
    scale_fill_manual(values = setNames(fills, fills), guide = "none") +
    scale_colour_manual(values = setNames(edges, edges), guide = "none") +
    # The axes are the per-run panels' in research/02_classify, which is what these
    # replace: still a log scale, but TICKED AND LABELLED AS THE EXPONENT, so a reader
    # reads "-2" rather than decoding 10^-2 or 1%. Decade ticks on y, every second decade
    # on x, and the titles name the transform.
    scale_x_log10(limits = XL, breaks = 10^seq(-6, 4, by = 2),
                  labels = seq(-6, 4, by = 2)) +
    scale_y_log10(limits = YL, breaks = 10^seq(-6, 0, by = 2),
                  labels = seq(-6, 0, by = 2)) +
    # NO AXIS TITLES. The four panels share one pair of axes in the composite and carry a
    # common title there, so a title on each would be set four times.
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.grid       = element_blank(),
          axis.line        = element_line(colour = "black", linewidth = 0.7),
          axis.ticks       = element_line(colour = "black", linewidth = 0.7),
          # the panels are placed small in the composite; the tick labels are the only
          # text left on them and have to survive that reduction
          axis.text        = element_text(colour = "black", size = 12),
          axis.title       = element_blank(),
          plot.margin      = margin(8, 10, 6, 6))
}

for (i in seq_len(nrow(pan))) {
  r <- pan$rung[i]
  f <- file.path(HERE, sprintf("runs_%d_%s.png", r,
                               gsub("_+$", "", gsub("[^a-z0-9]+", "_",
                                                    tolower(pan$label[i])))))
  ragg::agg_png(f, width = 3.2, height = 3.0, units = "in", res = 600,
                background = "white")
  print(panel_plot(r))
  invisible(dev.off())
  cat(sprintf("  %-34s %4d points, %2d called\n", basename(f),
              sum(pts$rung == r), sum(pts$called[pts$rung == r])))
}

# ── legend ────────────────────────────────────────────────────────────────────
# Keyed on the matrix SLOT, which is what the hue encodes. Keying it on taxon_class
# duplicated the key, because `resident` spans two slots and therefore two hues.
# The five categories come from runs_cats.tsv, which prep writes IN SCHEME ORDER and
# writes in full: a category with no point in these four runs is still part of the scheme
# the composite shares with samples.png, and dropping it would make the two legends differ
# for a reason the reader cannot see.
cats <- rd("runs_cats.tsv")
kd <- data.frame(k = factor(cats$category, levels = cats$category), x = 1, y = 1)
g1 <- ggplot(kd, aes(x, y, fill = k)) + geom_point(shape = 21, size = 2.6) +
  scale_fill_manual(values = setNames(cats$fill, cats$category), name = NULL,
                    guide = guide_legend(nrow = 1))

# The null is no longer drawn, so it is no longer a key. What remains is the two floors,
# which are the only lines on the panel and the only ones that decide a call.
lnames <- c(sprintf("abundance floor: %.2f x the null", SF),
            sprintf("coverage floor: %.0f%%", CF * 100))
ln <- data.frame(k = factor(lnames, levels = lnames), x = rep(1:2, each = 2),
                 y = rep(1:2, 2))
g2 <- ggplot(ln, aes(x, y, linetype = k, colour = k, group = k)) + geom_line() +
  scale_linetype_manual(values = c("22", "12"), name = NULL) +
  scale_colour_manual(values = c("#4a4a4a", "#4a4a4a"), name = NULL) +
  guides(linetype = guide_legend(nrow = 1), colour = guide_legend(nrow = 1))

# ONE ROW, as the samples legend is, so the two sit under the composite as a single band.
pick <- function(p) cowplot::get_plot_component(
  p + theme_void() + theme(legend.position = "bottom",
                           legend.text = element_text(size = 9.5, colour = "#1a1a1a",
                                                      margin = margin(l = 4, r = 14)),
                           legend.key.spacing.x = unit(2, "pt"),
                           legend.margin = margin(0, 0, 0, 6)),
  "guide-box-bottom", return_all = TRUE)
lg <- cowplot::plot_grid(pick(g1), pick(g2), nrow = 1, rel_widths = c(2.45, 1.55))
ragg::agg_png(file.path(HERE, "runs_legend.png"), width = 11.4, height = 0.32,
              units = "in", res = 600, background = "white")
print(lg)
invisible(dev.off())
cat(sprintf("wrote %d panels and runs_legend.png\n", nrow(pan)))
