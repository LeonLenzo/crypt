#!/usr/bin/env Rscript
# plot_runs.R — four columns of samples.png, opened up on the axes the call is made on.
#
#   python release/03_classify/figures/runs/prep_runs.py   # on Setonix; decides
#   Rscript release/03_classify/figures/runs/plot_runs.R   # draws
#
# Reads  runs_points.tsv, runs_floors.tsv, ../../data/inputs/panel_runs.tsv
# Writes runs.png, runs_legend.png
#
# Same split as the samples panel: prep decides, this draws. Every colour arrives as a hex
# computed from 02_reference's MATRIX, so a category is one colour across the whole chapter.
#
# ## The axes
#
# Redundancy available against coverage observed, both normalised by the taxon's own db2
# minimizer total, which is what makes the null ONE curve rather than one per taxon:
#
#   solid  y = 1-exp(-x)            what random sampling at that depth would give
#   dashed y = 0.05*(1-exp(-x))     the depth-normalised abundance floor
#   dotted y = 0.02                 the coverage floor, anchored on the off-host null
#
# Drawing the saturation floor as a CURVE rather than a horizontal line is the point of the
# 2026-11-24 measurement: it is not a shape test, it is an abundance floor that bends down
# at low depth, and the two floors meet at r = 0.51. Above that the dashed curve binds and
# below it the dotted line does, which the panel shows directly.

suppressPackageStartupMessages({library(ggplot2); library(grid); library(ggrepel)})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(HERE) || HERE == "") HERE <- "release/03_classify/figures/runs"
rd <- function(f, ...) read.delim(file.path(HERE, f), sep = "\t", quote = "",
                                  check.names = FALSE, colClasses = "character", ...)

pts <- rd("runs_points.tsv")
fl  <- rd("runs_floors.tsv")
pan <- read.delim(file.path(HERE, "../../data/inputs/panel_runs.tsv"), sep = "\t",
                  comment.char = "#", colClasses = "character")

num <- function(x) as.numeric(x)
pts$redundancy <- num(pts$redundancy); pts$coverage <- num(pts$coverage)
pts$rung <- as.integer(pts$rung)
gv <- function(k) num(fl$value[fl$item == k])
CF <- gv("coverage_floor"); SF <- gv("saturation_floor")

null <- fl[grepl("^null:", fl$item), ]
null <- data.frame(x = num(sub("^null:", "", null$item)), y = num(null$value))

pan$rung <- as.integer(pan$rung)
pan <- pan[order(pan$rung), ]
wrap <- function(x, n) vapply(strwrap(x, n, simplify = FALSE),
                              paste, "", collapse = "\n")
lab <- setNames(sprintf("%d  %s\n%s", pan$rung, pan$label, wrap(pan$expect, 46)), pan$rung)
pts$panel <- factor(lab[as.character(pts$rung)], levels = lab)

# A CALLED point is filled and ringed; a rejected one is hollow and pale. The rejected cloud
# is most of the data and carries no identity, so it is drawn small and unlabelled: what it
# shows is that the floors refuse far more than they admit, which is a shape, not a list.
pts$called <- pts$call == "detected"
fills <- sort(unique(pts$fill)); edges <- sort(unique(pts$edge))

p <- ggplot(pts, aes(redundancy, coverage)) +
  # no facet column on `null`, so ggplot repeats the curves in every panel. Assigning one
  # drew them in the first panel only, which made three of the four unreadable.
  geom_line(data = null, aes(x, y), inherit.aes = FALSE,
            colour = "#9a9a9a", linewidth = 0.5) +
  geom_line(data = null, aes(x, y * SF), inherit.aes = FALSE,
            colour = "#4a4a4a", linewidth = 0.5, linetype = "22") +
  geom_hline(yintercept = CF, colour = "#4a4a4a", linewidth = 0.4, linetype = "12") +
  geom_point(data = subset(pts, !called), shape = 21, size = 0.75, stroke = 0.2,
             fill = "#ffffff", colour = "#c4c4c4") +
  geom_point(data = subset(pts, called), aes(fill = fill, colour = edge),
             shape = 21, size = 2.1, stroke = 0.5) +
  # The field panel puts eleven labels in one corner, so repel needs room and force, and
  # the seed is fixed because an unseeded layout moves between renders.
  ggrepel::geom_text_repel(data = subset(pts, called), aes(label = taxon),
                           size = 2.05, fontface = "italic", colour = "#1a1a1a",
                           segment.colour = "#b0b0b0", segment.size = 0.25,
                           min.segment.length = 0, box.padding = 0.45,
                           point.padding = 0.15, force = 6, force_pull = 0.4,
                           max.overlaps = Inf, max.iter = 40000, max.time = 3,
                           seed = 1) +
  scale_fill_manual(values = setNames(fills, fills), guide = "none") +
  scale_colour_manual(values = setNames(edges, edges), guide = "none") +
  scale_x_log10(labels = function(x) format(x, scientific = TRUE, trim = TRUE)) +
  # a fixed-accuracy percent prints 0.0001 as "0.00%"; two significant figures keeps every
  # decade legible across the six the axis spans
  scale_y_log10(labels = function(y) paste0(signif(y * 100, 2), "%")) +
  facet_wrap(~ panel, nrow = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "redundancy available   (reads x k-mers per spot / the taxon's db2 minimizers)",
       y = "coverage observed   (distinct minimizers / the taxon's db2 minimizers)") +
  theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#f0f0f0", linewidth = 0.3),
        panel.border     = element_rect(fill = NA, colour = "#8f8f8f", linewidth = 0.45),
        panel.spacing    = unit(8, "pt"),
        strip.text       = element_text(size = 7.6, colour = "#1a1a1a", hjust = 0,
                                        lineheight = 1.25, margin = margin(b = 4)),
        axis.text        = element_text(size = 6.8, colour = "#707070"),
        axis.title       = element_text(size = 7.6, colour = "#707070"),
        plot.margin      = margin(6, 14, 6, 8))

ragg::agg_png(file.path(HERE, "runs.png"), width = 12.0, height = 4.6,
              units = "in", res = 600, background = "white")
print(p)
invisible(dev.off())
# ── legend ────────────────────────────────────────────────────────────────────
# Two keys: the point's reference-matrix category, which is the same hue as every other
# figure in the chapter, and the three drawn lines. Both are read from the data.
cats <- unique(pts[pts$called, c("taxon_class", "fill")])
NAMES <- c(wheat_pathogen = "wheat pathogen", plant_pathogen = "pathogen of another plant",
           resident = "resident", host = "host")
cats$nm <- NAMES[cats$taxon_class]
cats <- cats[order(match(cats$taxon_class, names(NAMES))), ]
kd <- data.frame(k = factor(cats$nm, levels = cats$nm), x = 1, y = 1)
g1 <- ggplot(kd, aes(x, y, fill = k)) + geom_point(shape = 21, size = 2.1) +
  scale_fill_manual(values = setNames(cats$fill, cats$nm), name = NULL,
                    guide = guide_legend(ncol = 1))

ln <- data.frame(k = factor(c("null: reads spread uniformly",
                              sprintf("abundance floor, %.2f x the null", SF),
                              sprintf("coverage floor, %.2f", CF)),
                            levels = c("null: reads spread uniformly",
                                       sprintf("abundance floor, %.2f x the null", SF),
                                       sprintf("coverage floor, %.2f", CF))),
                 x = 1:3, y = 1)
g2 <- ggplot(ln, aes(x, y, linetype = k, colour = k)) + geom_line() +
  scale_linetype_manual(values = c("solid", "22", "12"), name = NULL) +
  scale_colour_manual(values = c("#9a9a9a", "#4a4a4a", "#4a4a4a"), name = NULL,
                      guide = guide_legend(ncol = 1))

pick <- function(p) cowplot::get_plot_component(
  p + theme_void() + theme(legend.position = "right",
                           legend.text = element_text(size = 7.4, colour = "#1a1a1a",
                                                      margin = margin(l = 3)),
                           legend.key.spacing.y = unit(1, "pt")),
  "guide-box-right", return_all = TRUE)
lg <- cowplot::plot_grid(pick(g1), pick(g2), nrow = 1, rel_widths = c(1, 1.25))
ragg::agg_png(file.path(HERE, "runs_legend.png"), width = 5.2, height = 0.75,
              units = "in", res = 600, background = "white")
print(lg)
invisible(dev.off())

cat(sprintf("wrote %s/runs.png\n  %d points, %d called, %d panels\n",
            HERE, nrow(pts), sum(pts$called), nlevels(pts$panel)))
