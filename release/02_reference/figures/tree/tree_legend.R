#!/usr/bin/env Rscript
# tree_legend.R — the key for panel B, drawn as its own file.
#
#   Rscript release/02_reference/figures/tree/tree_legend.R
#
# Reads  tree_tips.tsv      (written by prep_tree.py, beside this script)
# Writes tree_legend.png
#
# Separate from the fan on purpose. A legend welded into the panel fixes its size and
# position relative to the plot, which makes the panel unplaceable; drawn apart, the key can
# sit wherever the layout wants.
#
# ONE SWATCH, ONE LABEL, one row each, following 01_cohort's cohort_network_legend. The
# inner ring's ramp is the sixth row, in the same grammar: its swatch is just a wider strip.
# Bar LENGTH is left to the caption, because the fan prints its own log gridline numerals
# and a key repeating them would be the only redundant thing here.
#
# The five colours are the five occupied observation cells of panel A. Keep in step with
# scripts/_palette.py MATRIX and figures/tree/tree.R CAT_COLOUR.

suppressPackageStartupMessages({ library(ggplot2) })

CAT <- c(wheat_pathogen   = "#e05a26",
         onwheat_pathogen = "#f5a800",
         onwheat_resident = "#1aa64b",
         plant_pathogen   = "#fbe49a",
         plant_resident   = "#c3e8bf")
# Wrapped to two lines where needed: the key is narrow because it shares the bottom band
# with the matrix, and a label that overruns is worse than one that takes a second line.
CAT_LABEL <- c(wheat_pathogen   = "wheat pathogen",
               onwheat_pathogen = "pathogen elsewhere,\nrecorded on wheat",
               onwheat_resident = "non-pathogen,\nrecorded on wheat",
               plant_pathogen   = "pathogen elsewhere,\nnot on wheat",
               plant_resident   = "non-pathogen,\nnot on wheat")

RING_PALE <- "#f2effa"; RING_DARK <- "#54278f"   # prep_tree.py ratio_colour() endpoints
OUT_IN <- 2.58                # fills the band beside the 4.89 in matrix: 7.67 - 4.89 - 0.2
SW <- 0.62                    # swatch side, in legend units
ROW <- 1.0
RAMP_W <- 2.3                 # the ring strip is a swatch, just a longer one
XMIN <- -0.1; XMAX <- 5.9     # extents are declared, so the PNG aspect can be derived
                              # from them rather than guessed; a stale guess clipped the
                              # ramp label off the right edge once already

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
tips <- read.delim(file.path(here, "tree_tips.tsv"), stringsAsFactors = FALSE)
stopifnot(all(unique(tips$category) %in% names(CAT)))   # a new category must reach the key

n <- length(CAT)
YMAX <- length(CAT) * 1.0 + 0.55
d <- data.frame(
  y    = seq(n, 1) * ROW,
  fill = unname(CAT),
  lab  = unname(CAT_LABEL))

# the inner-ring ramp, one row below the categories
ry <- -0.75
YMIN <- ry - 1.05
steps <- 60
edges <- seq(0, RAMP_W, length.out = steps + 1)
ramp <- data.frame(xmin = edges[-(steps + 1)], xmax = edges[-1],
                   fill = colorRampPalette(c(RING_DARK, RING_PALE))(steps))

p <- ggplot(d) +
  geom_rect(aes(xmin = 0, xmax = SW, ymin = y - SW / 2, ymax = y + SW / 2),
            fill = d$fill, colour = "#4a4a4a", linewidth = 0.45) +
  geom_text(aes(x = SW + 0.32, y = y, label = lab), hjust = 0, size = 3.1,
            colour = "#1a1a1a", lineheight = 0.92) +
  geom_rect(data = ramp, aes(xmin = xmin, xmax = xmax, ymin = ry - SW / 2,
                             ymax = ry + SW / 2), fill = ramp$fill, inherit.aes = FALSE) +
  annotate("rect", xmin = 0, xmax = RAMP_W, ymin = ry - SW / 2, ymax = ry + SW / 2,
           fill = NA, colour = "#4a4a4a", linewidth = 0.45) +
  annotate("text", x = 0, y = ry - SW / 2 - 0.28, label = "few", hjust = 0, size = 2.9,
           colour = "#707070") +
  annotate("text", x = RAMP_W, y = ry - SW / 2 - 0.28, label = "all", hjust = 1,
           size = 2.9, colour = "#707070") +
  annotate("text", x = 0, y = ry + SW / 2 + 0.34, label = "assemblies annotated",
           hjust = 0, size = 3.1, colour = "#1a1a1a") +
  coord_fixed(xlim = c(XMIN, XMAX), ylim = c(YMIN, YMAX), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(4, 4, 4, 4))

# coord_fixed means the PNG must carry the extents' own aspect or ggplot pads it with
# whitespace on one axis; derive it rather than hard-coding a ratio.
out <- file.path(here, "tree_legend.png")
ragg::agg_png(out, width = OUT_IN, height = OUT_IN * (YMAX - YMIN) / (XMAX - XMIN),
              units = "in", res = 300, background = "white")
print(p); invisible(dev.off())
cat(sprintf("  -> %s\n", out))
