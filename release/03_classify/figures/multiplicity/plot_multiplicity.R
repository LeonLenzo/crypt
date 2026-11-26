#!/usr/bin/env Rscript
# plot_multiplicity.R — draw the per-project panel from what prep_multiplicity.py decided.
#
#   python release/03_classify/figures/multiplicity/prep_multiplicity.py   # decides
#   Rscript release/03_classify/figures/multiplicity/plot_multiplicity.R   # draws
#
# Reads  multiplicity_bars.tsv      one row per bar segment: x range, colour
#        multiplicity_projects.tsv  one row per project: n, rate, depth, marker
#        multiplicity_ramp.tsv      the bins and the marker, for the legend
# Writes multiplicity.png, multiplicity_legend.png
#
# THIS SCRIPT DECIDES NOTHING. Each segment arrives as an x range and a hex string, so
# nothing here stacks a bar or picks a colour: a factor-level order cannot silently reverse
# a bar and the ramp cannot drift from the one the raster's legend shows.
#
# ## Two columns, one row per project
#
# Column one is the distribution: what share of the project's runs carry 0, 1, 2, 3, 4 or
# more wheat pathogens, as one bar. Column two is the median fungal reads per run, which is
# the depth that matters once the host is removed, and it is there to be checked against the
# first: the obvious reading of a low co-infection rate is a shallow library, and the two
# columns show that the projects with no co-infection are not the shallow ones.
#
# THE TWO COLUMNS SHARE ONE LINEAR x TRANSFORM and differ only in range, which is why the
# depth arrives pre-logged from prep_multiplicity.py. scale_x_log10 is a scale, not a panel
# setting, so it would have taken the share column with it.
#
# ## House style
#
# leon-figure-style: clean lines, little text, white background set explicitly. No axis
# titles -- the column strips say what each axis is -- and no panel borders: a bar has an
# edge of its own and a box around it is a second line saying the same thing. Rendered
# through ragg::agg_png for the same reason plot_samples.R is.

suppressPackageStartupMessages({library(ggplot2); library(grid); library(ggh4x)})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(HERE) || HERE == "") HERE <- "release/03_classify/figures/multiplicity"
rd <- function(f) read.delim(file.path(HERE, f), sep = "\t", quote = "",
                             check.names = FALSE, colClasses = "character")

bars <- rd("multiplicity_bars.tsv")
proj <- rd("multiplicity_projects.tsv")
ramp <- rd("multiplicity_ramp.tsv")
for (c in c("row", "n", "max", "virus")) proj[[c]] <- as.integer(proj[[c]])
for (c in c("rate", "mean", "log_f")) proj[[c]] <- as.numeric(proj[[c]])
for (c in c("x0", "x1", "share")) bars[[c]] <- as.numeric(bars[[c]])
bars$row <- as.integer(bars$row)

# Projects keep prep's order: block, then co-infection rate, then n. A discrete y scale runs
# bottom to top, so the levels are reversed and the highest rate sits at the top.
lv <- proj$project[order(proj$row)]
# The accessions get a PANEL of their own, the first facet column, the same device
# plot_samples.R uses for its row labels and for the same reason: drawn in the plot margin
# with clipping off they ran over the block strips, which sit between the margin and the
# first panel. A panel cannot collide with a strip.
LABEL <- "   "                    # a blank strip label, not a word
BARS  <- "wheat pathogens per run"
DEPTH <- "median fungal reads (M)"
fac <- function(d) {
  d$project <- factor(d$project, levels = rev(lv))
  d$block <- factor(d$block, levels = unique(proj$block[order(proj$row)]))
  d
}
bars <- fac(bars); proj <- fac(proj)
PARTS <- c(LABEL, BARS, DEPTH)
bars$part <- factor(BARS, levels = PARTS)
proj$part <- factor(DEPTH, levels = PARTS)

# The n goes past the right edge of the depth column with clipping off, the same idiom the
# raster uses for its row totals: a number in the margin is an annotation, not a panel.
# A CONSTANT SPACE PREFIX with hjust = 0 does the gap, because nudge_x would grow the
# panel's own scale and put the padding inside it.
DEPTH_HI <- 7.55
tot <- transform(proj, x = DEPTH_HI,
                 text = paste0(strrep(" ", 3), trimws(format(proj$n, big.mark = ","))))

# A VIRUS SURVEY'S ACCESSION IS DRAWN IN THE DECLARED GREY, not marked with a glyph of its
# own. Five of the six projects with no co-infection at all are virus surveys that declared
# no fungus because they were not looking for one, and the raster says that in the same
# grey, so the two panels make the point in one language.
#
# THE ACCESSIONS ARE A geom_text LAYER, NOT THE y AXIS. element_text takes a vector of
# colours only unofficially, and with two facet rows it recycles the vector per panel: the
# grey landed on five projects in the upper block that are not virus surveys at all. A text
# layer carries its colour with its row and cannot slip.
MARK <- ramp$fill[ramp$k == "virus"]
labs_df <- transform(proj, x = 1, part = factor(LABEL, levels = PARTS),
                     text = as.character(proj$project))
# The label panel needs an explicit 0..1 range, pinned by two invisible points, or its
# extent is wherever the text happens to sit and the column moves when a label changes.
pins <- data.frame(x = c(0, 1), project = proj$project[1], block = proj$block[1],
                   part = factor(LABEL, levels = PARTS))

fills <- sort(unique(bars$fill))

p <- ggplot(bars, aes(y = project)) +
  geom_tile(aes(x = (x0 + x1) / 2, width = x1 - x0, fill = fill), height = 0.74) +
  scale_fill_manual(values = setNames(fills, fills), guide = "none") +
  # The depth column: one point per project on the same rows. Hollow, so it reads as a
  # different kind of quantity from the filled bars beside it.
  geom_point(data = proj, aes(x = log_f), inherit.aes = TRUE, shape = 21,
             fill = "white", colour = "#1a1a1a", stroke = 0.6, size = 1.9) +
  geom_text(data = tot, aes(x = x, y = project, label = text), inherit.aes = FALSE,
            hjust = 0, size = 3.2, colour = "#3a3a3a") +
  geom_blank(data = pins, aes(x = x, y = project), inherit.aes = FALSE) +
  geom_text(data = labs_df, aes(x = x, y = project, label = text, colour = label_colour),
            inherit.aes = FALSE, hjust = 1, size = 3.3) +
  scale_colour_identity() +
  facet_grid(block ~ part, scales = "free", space = "free_y", switch = "y") +
  facetted_pos_scales(x = list(
    part == LABEL ~ scale_x_continuous(limits = c(0, 1), breaks = NULL,
                                       expand = expansion(mult = 0.01)),
    part == BARS ~ scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                                      labels = c("0", "", "50%", "", "100%"),
                                      expand = expansion(mult = 0.01)),
    # Breaks are decades of reads labelled in millions: the range spans 14.7 thousand to
    # 20.3 million, which is three orders of magnitude and cannot be read linearly.
    part == DEPTH ~ scale_x_continuous(limits = c(4, DEPTH_HI), breaks = 4:7,
                                       labels = c("0.01", "0.1", "1", "10"),
                                       expand = expansion(mult = 0.02)))) +
  force_panelsizes(cols = c(0.13, 0.65, 0.22)) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    # ONE grid line family: faint verticals in the panels, so a bar's end can be read
    # against 50% and a point against its decade without a box around either.
    panel.grid.major.x = element_line(colour = "#e4e4e4", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x  = unit(20, "pt"),
    panel.spacing.y  = unit(7, "pt"),
    panel.border     = element_blank(),
    axis.text.y      = element_blank(),
    axis.text.x      = element_text(size = 9, colour = "#3a3a3a"),
    axis.ticks       = element_blank(),
    strip.text.x     = element_text(size = 10.5, face = "bold", colour = "#1a1a1a",
                                    margin = margin(b = 5)),
    strip.text.y.left = element_text(size = 9.5, face = "bold", colour = "#3a3a3a",
                                     angle = 90, margin = margin(r = 4, l = 2)),
    strip.placement  = "outside",
    plot.margin      = margin(6, 34, 6, 6))

# Height is set from the content: one row per project plus the strips, so adding a project
# cannot squeeze the rows a draw later.
H <- 0.62 + 0.148 * nrow(proj) + 0.17 * nlevels(proj$block)
ragg::agg_png(file.path(HERE, "multiplicity.png"), width = 11.4, height = H,
              units = "in", res = 600, background = "white")
print(p)
invisible(dev.off())

# ── legend ────────────────────────────────────────────────────────────────────
# Two keys: the five bins, and what a grey accession means. Both are read from
# multiplicity_ramp.tsv, so a bin that stops occurring stops being drawn and nothing here
# can name a colour the panel did not use.
bins <- ramp[ramp$k != "virus", ]
kd <- data.frame(k = factor(bins$klabel, levels = bins$klabel), x = 1, y = 1)
g1 <- ggplot(kd, aes(x, y, fill = k)) + geom_tile() +
  scale_fill_manual(values = setNames(bins$fill, bins$klabel),
                    name = "wheat pathogens in the run  ",
                    guide = guide_legend(nrow = 1, keywidth = unit(11, "pt"),
                                         keyheight = unit(11, "pt")))
pick <- function(p) cowplot::get_plot_component(
  p + theme_void() + theme(legend.position = "bottom",
                           legend.title = element_text(size = 9.5, colour = "#1a1a1a"),
                           legend.text = element_text(size = 9.5, colour = "#1a1a1a",
                                                      margin = margin(l = 4, r = 14)),
                           legend.key.spacing.x = unit(2, "pt"),
                           legend.margin = margin(0, 0, 0, 6)),
  "guide-box-bottom", return_all = TRUE)

# THE VIRUS-SURVEY KEY IS A LINE OF GREY TEXT, not a swatch. The panel marks those projects
# by drawing their accession in the declared grey, and a key showing a grey square would be
# describing a mark the panel does not make.
note <- cowplot::ggdraw() +
  cowplot::draw_label(paste0("accessions in grey: ",
                             ramp$klabel[ramp$k == "virus"]),
                      colour = MARK, size = 9.5, x = 0, hjust = 0)

lg <- cowplot::plot_grid(pick(g1), note, nrow = 1, rel_widths = c(1.5, 1))
ragg::agg_png(file.path(HERE, "multiplicity_legend.png"), width = 11.4, height = 0.32,
              units = "in", res = 600, background = "white")
print(lg)
invisible(dev.off())

cat(sprintf("wrote %s/multiplicity.png  %.2f x %.2f in\n", HERE, 11.4, H))
cat(sprintf("  %d projects over %s runs, %d bar segments, rate %.1f%% to %.1f%%\n",
            nrow(proj), format(sum(proj$n), big.mark = ","), nrow(bars),
            100 * min(proj$rate), 100 * max(proj$rate)))
