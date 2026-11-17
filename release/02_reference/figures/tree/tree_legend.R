#!/usr/bin/env Rscript
# tree_legend.R — the key for panel B, drawn as its own file.
#
#   Rscript release/02_reference/figures/tree/tree_legend.R
#
# Reads  tree_tips.tsv, tree_arcs.tsv   (prep_tree.py and tree.R, beside this script)
# Writes tree_legend.png
#
# Separate from the fan on purpose. A legend welded into the panel fixes its size and
# position relative to the plot, which makes the panel unplaceable.
#
# TWO COLUMNS: bar colour on the left, arc symbol on the right, the ring ramp beneath them.
# One column would be ~5 in tall and the page band is 3.3; side by side it is 2.5.
#
# The SHAPE ASSIGNMENT IS READ BACK from tree_arcs.tsv, which tree.R writes while drawing,
# so the key cannot drift from the panel. The colours are mirrored by hand from
# scripts/_palette.py MATRIX and tree.R CAT_COLOUR; keep those three in step.

suppressPackageStartupMessages({ library(ggplot2) })

CAT <- c(wheat_pathogen   = "#e05a26",
         onwheat_pathogen = "#f5a800",
         onwheat_resident = "#1aa64b",
         plant_pathogen   = "#fbe49a",
         plant_resident   = "#c3e8bf")
CAT_LABEL <- c(wheat_pathogen   = "wheat pathogen",
               onwheat_pathogen = "pathogen elsewhere,\nrecorded on wheat",
               onwheat_resident = "non-pathogen,\nrecorded on wheat",
               plant_pathogen   = "pathogen elsewhere,\nnot on wheat",
               plant_resident   = "non-pathogen,\nnot on wheat")
RING_PALE <- "#f2effa"; RING_DARK <- "#54278f"   # prep_tree.py ratio_colour() endpoints
INK <- "#1a1a1a"; MUTED <- "#707070"

OUT_IN <- 3.67                # band is 7.67 wide: 3.80 matrix + 0.20 gutter + this
SW <- 0.58                    # swatch side, in legend units
ROW <- 1.05
C1 <- 0.0                     # column origins
C2 <- 6.35
RAMP_W <- 2.4
XMIN <- -0.15; XMAX <- 12.1   # declared, so the PNG aspect is derived and never clipped

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
tips <- read.delim(file.path(here, "tree_tips.tsv"), stringsAsFactors = FALSE)
arcs <- read.delim(file.path(here, "tree_arcs.tsv"), stringsAsFactors = FALSE)
stopifnot(all(unique(tips$category) %in% names(CAT)))   # a new category must reach the key

n <- length(CAT)
ytop <- n * ROW
col1 <- data.frame(y = seq(n, 1) * ROW, fill = unname(CAT), lab = unname(CAT_LABEL))
col2 <- data.frame(y = seq(nrow(arcs), 1) * ROW, shape = arcs$shape, lab = arcs$class)

ry <- -1.30                                     # the ramp row, under both columns
YMAX <- ytop + 1.05
YMIN <- ry - 1.00

ramp_edges <- seq(C1, C1 + RAMP_W, length.out = 61)
ramp <- data.frame(xmin = ramp_edges[-61], xmax = ramp_edges[-1],
                   fill = colorRampPalette(c(RING_DARK, RING_PALE))(60))

p <- ggplot() +
  # headers
  annotate("text", x = C1, y = ytop + 0.62, label = "bar colour", hjust = 0, size = 3.4,
           fontface = "bold", colour = INK) +
  annotate("text", x = C2, y = ytop + 0.62, label = "arc symbol", hjust = 0, size = 3.4,
           fontface = "bold", colour = INK) +
  # column 1: the five panel-A cells
  geom_rect(data = col1, aes(xmin = C1, xmax = C1 + SW, ymin = y - SW / 2, ymax = y + SW / 2),
            fill = col1$fill, colour = "#4a4a4a", linewidth = 0.45) +
  geom_text(data = col1, aes(x = C1 + SW + 0.3, y = y, label = lab), hjust = 0, size = 3.0,
            colour = INK, lineheight = 0.92) +
  # column 2: the taxonomic classes marked on the fan
  geom_point(data = col2, aes(x = C2 + SW / 2, y = y, shape = shape), size = 2.9,
             fill = "white", colour = "#333333", stroke = 0.75) +
  scale_shape_identity() +
  geom_text(data = col2, aes(x = C2 + SW + 0.3, y = y, label = lab), hjust = 0, size = 3.0,
            colour = INK) +
  # the annotation-ratio ramp
  annotate("text", x = C1, y = ry + SW / 2 + 0.40, label = "assemblies annotated (inner ring)",
           hjust = 0, size = 3.0, colour = INK) +
  geom_rect(data = ramp, aes(xmin = xmin, xmax = xmax, ymin = ry - SW / 2, ymax = ry + SW / 2),
            fill = ramp$fill) +
  annotate("rect", xmin = C1, xmax = C1 + RAMP_W, ymin = ry - SW / 2, ymax = ry + SW / 2,
           fill = NA, colour = "#4a4a4a", linewidth = 0.45) +
  annotate("text", x = C1, y = ry - SW / 2 - 0.30, label = "few", hjust = 0, size = 2.8,
           colour = MUTED) +
  annotate("text", x = C1 + RAMP_W, y = ry - SW / 2 - 0.30, label = "all", hjust = 1,
           size = 2.8, colour = MUTED) +
  coord_fixed(xlim = c(XMIN, XMAX), ylim = c(YMIN, YMAX), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(4, 4, 4, 4))

# coord_fixed means the PNG must carry the extents' own aspect, or ggplot pads one axis.
out <- file.path(here, "tree_legend.png")
ragg::agg_png(out, width = OUT_IN, height = OUT_IN * (YMAX - YMIN) / (XMAX - XMIN),
              units = "in", res = 300, background = "white")
print(p); invisible(dev.off())
cat(sprintf("  -> %s  (%d colours, %d classes)\n", out, nrow(col1), nrow(col2)))
