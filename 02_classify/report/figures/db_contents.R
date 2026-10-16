#!/usr/bin/env Rscript
# db_contents.R — what db_v4 can see, as a circular taxonomic cladogram.
#
#   Rscript 02_classify/report/figures/db_contents.R                     # assembly mode
#   Rscript 02_classify/report/figures/db_contents.R --self-contained    # + title, legend
#
# Reads only what prep_db_contents.py wrote. Run that first; it owns the membership logic,
# the taxonomy walk and the rank ladder, so this file is layout and style and nothing else.
#
# IT IS A TAXONOMY, NOT A PHYLOGENY. Branch lengths are all 1 because every edge is one rank
# step, so radius encodes taxonomic RANK (kingdom -> species) and not divergence. Never
# caption it as a phylogeny.
#
# Why ggtree and not ggtreeExtra: the fan layout is coord_polar(theta = "y") under the hood,
# so a plain geom_rect in tree coordinates becomes a WEDGE that tiles its angular slot
# exactly. That is the property the bars need -- a fixed-width line would be wider than its
# slot at this tip count and would overdraw its neighbours -- and it comes for free.

suppressPackageStartupMessages({
  library(ggplot2); library(ggtree); library(treeio); library(ape); library(dplyr)
})

# ── everything tunable ────────────────────────────────────────────────────────
ROUTE_COLOUR <- c(
  host            = "#3e9e3b",
  target_cereal   = "#8c0404",
  target_other    = "#e77e7e",
  host_associate  = "#22127e",
  plant_associate = "#add0e4",
  congener        = "#dfd8d8")
ROUTE_LABEL <- c(
  target_cereal   = "pathogen recorded on wheat, rice, maize or sorghum",
  target_other    = "plant pathogen, not recorded on those crops",
  host_associate  = "host associate: seen in a phytobiome study",
  plant_associate = "plant associate: recorded on a plant",
  congener        = "congener: shares a genus with a pathogen",
  host            = "host (plant)")

BRANCH_COLOUR     <- "#8a8a8a"   # internal edges
BRANCH_WIDTH      <- 0.12
TIP_BRANCH_COLOUR <- "#8a8a8a"   # terminal edges; see the note below
TIP_BRANCH_WIDTH  <- 0.06

HOLE       <- 5    # tree depth units of empty centre. BIGGER = bigger hole AND shorter
                     # branches, because the 7 rank rings are squeezed into what is left.
BAR_GAP    <- 0.35   # gap between the tip ring and the start of the bars
BAR_SCALE  <- 2   # radial length of one log10 decade, in tree depth units
ARC_GAP    <- 0.45   # gap between the longest bar and the clade arcs
ARC_WIDTH  <- 1.1    # weight of the clade arcs
ARC_PAD    <- 5      # tips trimmed from EACH END of a clade arc, so neighbours separate.
                     # Without it the arcs run tip to tip and consecutive classes sit one tip
                     # apart -- 0.22 degrees at this tip count, which reads as one continuous
                     # ring. Capped at 15% of the clade so a small one keeps a visible arc.
LABEL_MIN  <- 40     # a class needs this many species to get an arc

# ── inputs ────────────────────────────────────────────────────────────────────
# base read.delim, NOT readr: vroom fails to load on this machine with
# "undefined symbol: SET_GROWABLE_BIT".
here  <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "02_classify/report/figures"
stem  <- file.path(here, "db_contents")
self_contained <- "--self-contained" %in% commandArgs(TRUE)

# The gap angle is NOT set here. prep_db_contents.py owns it, because the clade angles it
# reports for Illustrator are computed from it; declaring it in both places meant a change on
# this side would silently leave the caption describing a different figure. To widen the gap,
# change GAP_DEG in prep_db_contents.py and re-run it.
params     <- read.delim(paste0(stem, "_params.tsv"), stringsAsFactors = FALSE)
OPEN_ANGLE <- params$value[params$param == "open_angle_deg"]
stopifnot(length(OPEN_ANGLE) == 1, is.finite(OPEN_ANGLE))

tr   <- read.tree(paste0(stem, "_tree.nwk"))
tips <- read.delim(paste0(stem, "_tips.tsv"), stringsAsFactors = FALSE)
tips$route <- factor(tips$route, levels = names(ROUTE_COLOUR))

cat(sprintf("  %d tips, %d internal nodes\n", Ntip(tr), Nnode(tr)))

# ── the tree ──────────────────────────────────────────────────────────────────
# Terminal and internal edges are drawn as two layers on purpose. There are as many terminal
# edges as species, 0.22 degrees apart, so at the internal weight they merge into a solid ring
# and bury the structure the figure exists to show.
# `color` and `linewidth`, NOT `colour` and `size`: ggtree forwards ... to the tree layer,
# and only that spelling reaches it. With `colour`/`size` the arguments are silently swallowed
# into fortify() and the branch settings do nothing -- the figure still draws, in defaults, so
# there is nothing to notice. ggtree warns about its own internal deprecations either way.
p <- suppressWarnings(
  ggtree(tr, layout = "fan", open.angle = OPEN_ANGLE,
         color = BRANCH_COLOUR, linewidth = BRANCH_WIDTH))

d       <- p$data
is_tip  <- d$isTip
tip_seg <- d[is_tip, ]
tip_seg$xpar <- d$x[match(tip_seg$parent, d$node)]
p <- p + geom_segment(data = tip_seg,
                      aes(x = xpar, xend = x, y = y, yend = y),
                      colour = TIP_BRANCH_COLOUR, linewidth = TIP_BRANCH_WIDTH,
                      inherit.aes = FALSE)

# ── the bars ──────────────────────────────────────────────────────────────────
xtip <- max(d$x[is_tip])
bars <- d[is_tip, c("label", "y")] |> left_join(tips, by = "label")
stopifnot(!anyNA(bars$route))
bars$x0 <- xtip + BAR_GAP
bars$x1 <- bars$x0 + log10(bars$n_annotated + 1) * BAR_SCALE

p <- p + geom_rect(data = bars,
                   aes(xmin = x0, xmax = x1, ymin = y - 0.5, ymax = y + 0.5, fill = route),
                   inherit.aes = FALSE) +
  scale_fill_manual(values = ROUTE_COLOUR, labels = ROUTE_LABEL,
                    breaks = names(ROUTE_LABEL), name = NULL, drop = FALSE)

# ── log gridlines, and the scale numerals that stay in BOTH modes ─────────────
# Under coord_polar(theta = "y") a constant x is a circle, so geom_vline draws the rings.
top   <- max(bars$n_annotated)
ticks <- c(1, 10, 100, 1000); ticks <- ticks[ticks <= top * 1.5]
tickx <- bars$x0[1] + log10(ticks + 1) * BAR_SCALE
p <- p + geom_vline(xintercept = tickx, colour = "#b0b0b0",
                    linewidth = 0.2, linetype = "dotted")

# The numerals sit in the open wedge, which is at the END of the y range.
ylab <- max(d$y) + Ntip(tr) * OPEN_ANGLE / 360 / 2
p <- p + annotate("text", x = tickx, y = ylab, label = format(ticks, big.mark = ","),
                  size = 4.2, colour = "#333333", hjust = 0.5, vjust = -0.4)

# ── clade arcs: marks only, never text. Positions are in db_contents_clades.tsv ──
arcs <- bars |> filter(class != "") |> group_by(class) |>
  summarise(n = n(), ymin = min(y), ymax = max(y), .groups = "drop") |>
  filter(n >= LABEL_MIN) |>
  mutate(pad  = pmin(ARC_PAD, 0.15 * n),
         ymin = ymin + pad, ymax = ymax - pad)
arcx <- bars$x0[1] + log10(top + 1) * BAR_SCALE + ARC_GAP
p <- p + geom_segment(data = arcs,
                      aes(x = arcx, xend = arcx, y = ymin, yend = ymax),
                      colour = "#333333", linewidth = ARC_WIDTH, lineend = "butt",
                      inherit.aes = FALSE)
cat(sprintf("  %d class arcs, smallest gap between neighbours %.1f deg\n",
            nrow(arcs), 2 * min(arcs$pad) * (360 - OPEN_ANGLE) / Ntip(tr)))

# ── the hole, the theme, and out ──────────────────────────────────────────────
# A NEGATIVE lower xlim pushes the root off the centre, which is how a fan layout gets its
# hole. Bigger HOLE also shortens every branch, since the rank rings share what is left.
p <- p + xlim(c(-HOLE, arcx + 0.1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        plot.margin      = margin(6, 6, 6, 6))

if (self_contained) {
  p <- p + labs(
    title    = "What db_v4 can see",
    subtitle = sprintf(paste0("%s fungal species carrying a gene model, so usable in a ",
                              "CDS-only index.  Bar: annotated assemblies, log scale\n",
                              "colour is the route into scope"),
                       format(nrow(bars), big.mark = ","))) +
    theme(legend.position = "bottom", legend.text = element_text(size = 11),
          plot.title = element_text(size = 20, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5, colour = "#333333")) +
    guides(fill = guide_legend(ncol = 2, override.aes = list(linewidth = 0)))
} else {
  p <- p + theme(legend.position = "none")
}

# ggtree opens the fan at 3 o'clock, so rotate it to the top where a reader looks first.
# The extra OPEN_ANGLE/2 is not a fudge: rotate_tree sets where the gap's leading EDGE lands,
# not its centre, so rotating by 90 alone leaves the wedge sitting half its own width
# clockwise of twelve. Applied LAST, because it rewrites the coord and any layer added
# afterwards would not follow it.
p <- rotate_tree(p, 90 + OPEN_ANGLE / 2)

png_path <- paste0(stem, ".png"); pdf_path <- paste0(stem, ".pdf")
ragg::agg_png(png_path, width = 15, height = if (self_contained) 16.2 else 15,
              units = "in", res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(pdf_path, width = 15, height = if (self_contained) 16.2 else 15, bg = "white")
print(p); invisible(dev.off())
cat(sprintf("  -> %s / .pdf   (mode: %s)\n", png_path,
            if (self_contained) "self-contained" else "assembly"))
