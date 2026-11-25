#!/usr/bin/env Rscript
# tree.R — the reference contents as a circular taxonomic cladogram (panel B).
#
#   Rscript release/02_reference/figures/tree/tree.R
#
# The panel carries NO legend and no title: tree_legend.R draws the key as its own file, so
# the fan can be placed at whatever size the layout wants without dragging a legend block
# with it. caption.py writes the words.
#
# Reads only what prep_tree.py wrote. IT IS A TAXONOMY, NOT A PHYLOGENY: radius encodes rank
# (kingdom -> species), not divergence. Bar length is the species' genome count on a log
# scale; tip colour is the label category. Layout and style only -- prep_tree.py owns the data.
#
# Structure follows 02_classify's db_contents.R: the fan is coord_polar(theta="y") under the
# hood, so a geom_rect in tree coordinates tiles its angular slot exactly (a wedge), which is
# what the bars need at this tip count.

suppressPackageStartupMessages({
  library(ggplot2); library(ggtree); library(treeio); library(ape); library(dplyr)
})

# ── style ─────────────────────────────────────────────────────────────────────
# FIVE colours, one per occupied observation cell of panel A, and the same five hexes, so a
# hue means the same thing in both panels. Residents are the recessive pale mass and the
# pathogens the few bright tips, which is the point: the eye reads the pathogen sliver
# against the community. Order here is the DRAW order, residents underneath.
# Keep in step with scripts/_palette.py MATRIX and figures/tree/tree_legend.R.
CAT_COLOUR <- c(
  wheat_pathogen   = "#e05a26",   # Wheat x Wheat pathogen
  onwheat_pathogen = "#f5a800",   # Wheat x Plant pathogen
  onwheat_resident = "#1aa64b",   # Wheat x Non pathogen
  plant_pathogen   = "#fbe49a",   # Plant x Plant pathogen
  plant_resident   = "#c3e8bf")   # Plant x Non pathogen

# OUT_IN is the size the panel is PLACED at, and every stroke and type size below is scaled
# by it. The design was tuned at 15 in; rendering at 15 and shrinking the PNG into an A4
# layout would take 12 pt numerals down to about 5 pt, so the panel is rendered at its final
# size instead and the constants scale with it. A4 landscape (11.69 x 8.27 in) holds this
# A4 PORTRAIT (8.27 x 11.69 in), fan on top and the matrix band beneath. The fan is limited
# by page WIDTH here, not by the matrix: 8.27 - 0.6 margins = 7.67, taken as 7.5 to leave the
# band its 3.26 + 0.2 gutter inside the 11.69 height.
OUT_IN <- 7.5
S <- OUT_IN / 15

BRANCH_COLOUR <- "#8a8a8a"; BRANCH_WIDTH <- 0.12 * S
TIP_BRANCH_COLOUR <- "#8a8a8a"; TIP_BRANCH_WIDTH <- 0.06 * S
HOLE <- 5; BAR_GAP <- 0.62; BAR_SCALE <- 2      # tree coordinates: NOT scaled by S
ARC_GAP <- 0.45; ARC_WIDTH <- 1.1 * S; ARC_PAD <- 5; LABEL_MIN <- 40
PRED_COLOUR <- "#1a1a1a"                   # the predicted-model ring, outside the band
SYM_GAP <- 1.05            # symbol centre, outside the arc. In TREE x-units, so not * S:
                           # the fan's radius is in the same units and both scale together.
# Filled shapes with a separate border, so they read as marks and not as data points. In
# class-name order; tree_legend.R reads the assignment back from tree_arcs.tsv.
SHAPES <- c(21, 22, 23, 24, 25)
TICK_SIZE <- 6.0 * S        # the bar-scale numerals; 4.2 printed at about 6 pt

# ── inputs ────────────────────────────────────────────────────────────────────
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
stem <- normalizePath(file.path(here, "tree"), mustWork = FALSE)
params <- read.delim(paste0(stem, "_params.tsv"), stringsAsFactors = FALSE)
OPEN_ANGLE <- params$value[params$param == "open_angle_deg"]
tr   <- read.tree(paste0(stem, "_tree.nwk"))
tips <- read.delim(paste0(stem, "_tips.tsv"), stringsAsFactors = FALSE)
tips$category <- factor(tips$category, levels = names(CAT_COLOUR))
cat(sprintf("  %d tips, %d internal nodes\n", Ntip(tr), Nnode(tr)))

# ── tree ──────────────────────────────────────────────────────────────────────
p <- suppressWarnings(ggtree(tr, layout = "fan", open.angle = OPEN_ANGLE,
                             color = BRANCH_COLOUR, linewidth = BRANCH_WIDTH))
d <- p$data; is_tip <- d$isTip
tip_seg <- d[is_tip, ]; tip_seg$xpar <- d$x[match(tip_seg$parent, d$node)]
p <- p + geom_segment(data = tip_seg, aes(x = xpar, xend = x, y = y, yend = y),
                      colour = TIP_BRANCH_COLOUR, linewidth = TIP_BRANCH_WIDTH,
                      inherit.aes = FALSE)

# ── bars: genome count, log scale, coloured by category ───────────────────────
xtip <- max(d$x[is_tip])
bars <- d[is_tip, c("label", "y", "angle")] |> left_join(tips, by = "label")
stopifnot(!anyNA(bars$category))
# draw residents first so the bright categories sit on top where arcs overlap
bars <- bars[order(bars$category, decreasing = TRUE), ]
bars$x0 <- xtip + BAR_GAP
bars$x1 <- bars$x0 + log10(bars$n_genomes + 1) * BAR_SCALE
p <- p + geom_rect(data = bars,
                   aes(xmin = x0, xmax = x1, ymin = y - 0.5, ymax = y + 0.5, fill = category),
                   inherit.aes = FALSE) +
  scale_fill_manual(values = CAT_COLOUR, breaks = names(CAT_COLOUR),
                    name = NULL, drop = FALSE)

# ── annotation-ratio band: a ring at the tip base under the bars ───────────────
# Per species, annotated assemblies / total assemblies. Pale = every assembly is annotated;
# dark purple = few are, i.e. genomes sitting unannotated that gene prediction would unlock.
# Drawn with a fixed fill VECTOR (not an aes) so it needs no scale and never collides with the
# category scale above. band_colour is precomputed in prep_tree.py.
band_x0 <- xtip + 0.04; band_x1 <- xtip + 0.30
p <- p + geom_rect(data = bars, aes(ymin = y - 0.5, ymax = y + 0.5),
                   xmin = band_x0, xmax = band_x1, fill = bars$band_colour,
                   inherit.aes = FALSE)

# ── predicted gene models: a second ring, just outside the annotation band ────
# One tick per tip whose gene model was PREDICTED here rather than downloaded from NCBI.
# Every one of those taxa had no NCBI annotation at all, so a tick means the whole bar
# beside it is ours and the species was uncallable before db2.
#
# It sits at the TIP BASE beside the band and not out past the bars, where a first attempt
# put it: beyond the longest bar is an annulus most tips never reach, so the marks floated
# unattached to anything and read as dust. Two thin rings at the tip base read as one small
# block per species, which is what they are: what NCBI annotated, then what we did. It is
# also kept out of the bar itself, because the bar carries the panel-A category hue and
# recolouring it would cost the figure its one colour vocabulary. Adjacent ticks merge into
# a solid arc on purpose, so a systematically unannotated clade reads as one.
pred <- bars[bars$n_predicted > 0, ]
pred_x0 <- xtip + 0.34; pred_x1 <- xtip + BAR_GAP - 0.04
p <- p + geom_rect(data = pred, aes(ymin = y - 0.5, ymax = y + 0.5),
                   xmin = pred_x0, xmax = pred_x1, fill = PRED_COLOUR,
                   inherit.aes = FALSE)
cat(sprintf("  %d tips carry a predicted gene model\n", nrow(pred)))

# ── log gridlines + numerals (kept in both modes) ─────────────────────────────
top <- max(bars$n_genomes)
ticks <- c(1, 10, 100, 1000); ticks <- ticks[ticks <= top * 1.5]
tickx <- bars$x0[1] + log10(ticks + 1) * BAR_SCALE
p <- p + geom_vline(xintercept = tickx, colour = "#b0b0b0", linewidth = 0.2 * S,
                    linetype = "dotted")
ylab <- max(d$y) + Ntip(tr) * OPEN_ANGLE / 360 / 2
p <- p + annotate("text", x = tickx, y = ylab, label = format(ticks, big.mark = ","),
                  size = TICK_SIZE, colour = "#333333", hjust = 0.5, vjust = -0.4)

# ── clade arcs, MARKED WITH A SYMBOL ──────────────────────────────────────────
# A symbol, not a name. Set radially, a 15-character class name is the least legible thing
# on the panel: it reads at a different angle in every quadrant, half of it needs flipping,
# and at this size it competes with the bars it is meant to annotate. A point shape is
# orientation-free, reads at a glance, and the name moves to the key, where there is room
# for it. The shape assignment is written to tree_arcs.tsv so tree_legend.R cannot drift
# from it.
ROT <- 90 + OPEN_ANGLE / 2
arcs <- bars |> filter(class != "") |> group_by(class) |>
  summarise(n = n(), ymin = min(y), ymax = max(y), .groups = "drop") |>
  filter(n >= LABEL_MIN) |>
  arrange(class) |>
  mutate(pad = pmin(ARC_PAD, 0.15 * n), ymin = ymin + pad, ymax = ymax - pad,
         ymid = (ymin + ymax) / 2,
         shape = SHAPES[seq_len(n())])
stopifnot(nrow(arcs) <= length(SHAPES))     # more classes than shapes = silent recycling
arcx <- bars$x0[1] + log10(top + 1) * BAR_SCALE + ARC_GAP
p <- p + geom_segment(data = arcs, aes(x = arcx, xend = arcx, y = ymin, yend = ymax),
                      colour = "#333333", linewidth = ARC_WIDTH, lineend = "butt",
                      inherit.aes = FALSE) +
  geom_point(data = arcs, aes(x = arcx + SYM_GAP, y = ymid, shape = shape),
             size = 7.5 * S, fill = "white", colour = "#333333", stroke = 1.6 * S,
             inherit.aes = FALSE) +
  scale_shape_identity()
write.table(arcs[, c("class", "n", "shape")], paste0(stem, "_arcs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("  %d class arcs: %s\n", nrow(arcs),
            paste(sprintf("%s(%d)", arcs$class, arcs$shape), collapse = ", ")))

# ── hole, theme, out ──────────────────────────────────────────────────────────
p <- p + xlim(c(-HOLE, arcx + SYM_GAP + 0.6)) +   # room for the symbols outside the arcs
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 0, 0, 0))

p <- p + theme(legend.position = "none")
p <- rotate_tree(p, ROT)

png <- paste0(stem, ".png")
ragg::agg_png(png, width = OUT_IN, height = OUT_IN, units = "in", res = 300,
              background = "white")
print(p); invisible(dev.off())
cat(sprintf("  -> %s\n", png))
