#!/usr/bin/env Rscript
# plot_samples.R — draw the per-run panel from what prep_samples.py decided.
#
#   python release/03_classify/figures/samples/prep_samples.py   # decides
#   Rscript release/03_classify/figures/samples/plot_samples.R   # draws
#
# Reads  samples_cells.tsv  one row per drawn cell: column, row, state, fill
#        samples_rows.tsv   one row per taxon row: label, band, half, colour
#        samples_cols.tsv   one row per run column: region, country, project
# Writes samples.png, samples_legend.png
#
# THIS SCRIPT DECIDES NOTHING. Every colour arrives as a hex string and every grouping as a
# factor level, both computed in prep_samples.py from the pipeline's own functions. If a
# question here needs an answer that is not already in the three tables, it belongs upstream.
#
# ## Why this half is R
#
# The panel is a grid of blocks: rows banded by disease or by reference slot, columns blocked
# by region, with gaps between both. In matplotlib that was hand-computed inches -- block
# extents, label collision tiers, rules positioned in gaps -- and most of the script was
# arithmetic about where things go. `facet_grid(space = "free")` is exactly that layout, so
# the arithmetic disappears and the panel widths stay proportional to the runs they hold.
#
# ## House style
#
# leon-figure-style: clean lines, little text, white background set explicitly, heavy strokes.
# Rendered through ragg::agg_png because the default device substitutes glyphs (the dagger on
# a flagged row is one). `read.delim`, not readr, which does not load on this machine
# (vroom.so, undefined symbol: SET_GROWABLE_BIT).

suppressPackageStartupMessages({library(ggplot2); library(grid); library(ggh4x)})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(HERE) || HERE == "") HERE <- "release/03_classify/figures/samples"
rd <- function(f) read.delim(file.path(HERE, f), sep = "\t", quote = "",
                             check.names = FALSE, colClasses = "character")

cells <- rd("samples_cells.tsv")
rows  <- rd("samples_rows.tsv")
cols  <- rd("samples_cols.tsv")
cells$col <- as.integer(cells$col); cells$row <- as.integer(cells$row)
rows$row  <- as.integer(rows$row);  cols$col  <- as.integer(cols$col)

# Regions keep the west-to-east order prep_samples.py put them in, which is the mean
# longitude of each region's own runs. Rows keep its order too: band, then prevalence.
cols$region <- factor(cols$region, levels = unique(cols$region[order(cols$col)]))
rows$band   <- factor(rows$band,   levels = unique(rows$band[order(rows$row)]))
rows$half   <- factor(rows$half,   levels = unique(rows$half[order(rows$row)]))
rows$label  <- factor(rows$label,  levels = rev(rows$label[order(rows$row)]))

d <- merge(merge(cells, rows[, c("row", "label", "band", "half")], by = "row"),
           cols[, c("col", "region")], by = "col")

# THE ROW LABELS GET THEIR OWN PANEL, so they sit in a bordered box like every section of
# the grid rather than floating in the axis margin. It is an extra facet column holding one
# text per row and no tiles, placed first, with a blank column strip.
#
# That means panel widths can no longer come from `space = "free"`, which sizes a panel by
# its data range and would make a one-unit-wide label panel invisible. ggh4x's
# force_panelsizes takes explicit relative widths instead: the label panel gets a fixed
# share and each region keeps its true share of the rest, so a region is still as wide as
# the runs it holds.
LABEL_PANEL <- "  "                    # a blank strip label, not a word
nreg <- table(cols$region)
d$region  <- factor(as.character(d$region), levels = c(LABEL_PANEL, names(nreg)))
# ROW TOTALS get a panel of their own on the right, for the same reason the labels do: a
# number in a bordered box is part of the grid, a number floating past the last region is
# not. n_runs comes from prep_samples.py, counted over the drawn columns.
lev <- c(LABEL_PANEL, names(nreg))
# The label panel gets an explicit 0..1 x range, pinned by two invisible points, so the two
# things inside it can be placed: the category at the left, the taxa in a column beside it.
# Without the pins the range is whatever the text happens to sit at and the layout moves
# when a label changes.
# The two columns moved apart when the text grew: at 3.6 pt "Blotches" ran into
# "P. nodorum". The split is set from the longest band name, not nudged until it looks
# right, so adding a band cannot silently re-create the collision.
LAB_X_CAT  <- 0.02                     # category name, left edge
LAB_X_TAXA <- 0.47                     # taxa column, left-justified
# A band name wider than the column runs under the binomials. Wrapping it is better than
# widening the column, which would cost grid width on every row to fit one word.
labs <- data.frame(x = LAB_X_TAXA, label = rows$label, band = rows$band, half = rows$half,
                   region = factor(LABEL_PANEL, levels = lev),
                   text = as.character(rows$label),
                   face = ifelse(rows$italic == "1", "italic", "bold"))
pins <- data.frame(x = c(0, 1), label = rows$label[1], band = rows$band[1],
                   half = rows$half[1], region = factor(LABEL_PANEL, levels = lev))
# The BAND NAME lives in the label box too, not in a strip beside it, so one box carries the
# category and the taxa it contains. One row per band, placed at the box's left edge and
# centred on the band's rows; the box is tinted with the band's prevailing matrix colour.
# The category is placed on the band's FIRST row, Leon's call. A discrete y scale has no
# "centre" to address -- each band panel holds only its own labels, so the old midpoint was
# the label at ceiling(k/2), which lands a half-row high on an even count and made the
# category drift up and down relative to the species beside it as bands changed size.
# The first row is exact, it is the same place in every band, and a left-aligned column of
# headings reading down the panel is easier to scan than five differently-offset ones.
#
# A WRAPPED name is anchored by its LINE COUNT, not pinned to row 1. geom_text centres a
# multi-line block on its y, so "Spots and\nBlotches" pinned to row 1 put its first line
# above the panel. Anchoring at ceiling(nlines / 2) lands line 1 on row 1 and line 2 on
# row 2, which is top-aligned for any number of lines without touching vjust.
# ONE geom_text PER WRAPPED LINE, each centred on its own species row. A multi-line label
# drawn as a single string is centred on the block, so its first line never sits level with
# the first species: vjust and nudge_y only approximate the offset, and the approximation
# moves with the font size. Splitting the name and placing line i on row i is exact, needs
# no nudge, and generalises to any number of lines.
bl <- rows[!duplicated(rows[, c("half", "band")]), ]
bands_df <- do.call(rbind, Map(function(h, b) {
  lb <- rows$label[rows$half == h & rows$band == b]
  # as.character: mapply/Map over a factor column hands the callback a one-element FACTOR,
  # and strwrap on that silently returns its integer code
  parts <- strwrap(as.character(b), width = 12)
  # A band with fewer rows than the name has lines cannot show the wrap -- "Virus survey"
  # over a one-row band lost the word "survey" entirely. Fall back to the unwrapped name
  # and let it be the widest thing in the column; a truncated label is the worse failure.
  if (length(parts) > length(lb)) parts <- as.character(b)
  n <- min(length(parts), length(lb))
  data.frame(band = factor(b, levels = levels(rows$band)),
             half = factor(h, levels = levels(rows$half)),
             region = factor(LABEL_PANEL, levels = lev),
             text = parts[seq_len(n)],
             ymid = factor(as.character(lb[seq_len(n)]), levels = levels(rows$label)))
}, bl$half, bl$band))

frames <- expand.grid(band = factor(levels(rows$band), levels = levels(rows$band)),
                      region = factor(setdiff(lev, LABEL_PANEL), levels = lev))
LAB_SHARE <- 0.26                      # of the drawn width, enough for the longest binomial
widths <- c(LAB_SHARE, as.numeric(nreg) / sum(as.numeric(nreg)) * (1 - LAB_SHARE))

# ROW TOTALS go on a secondary discrete axis, NOT in a panel of their own. A panel would
# carry the grid's border and these are an annotation, not a section of the heatmap.
# ROW TOTALS are drawn past the right edge of the last region panel with clipping off,
# NOT as a panel (a panel carries the grid's border and these are an annotation) and NOT on
# a secondary axis (with free scales the duplicated axis resolved its labels against the
# wrong panel and printed NA).
# The gap before the totals is a CONSTANT SPACE PREFIX with hjust = 0, which is the only
# version of this that both left-aligns and leaves the scale alone. Two others were tried:
# nudge_x is part of the positioned data, so the panel's x scale grew to contain it and the
# padding appeared INSIDE the box; a negative hjust offsets by a fraction of each string's
# OWN width, so "1,029" and "96" started in different places.
last_reg <- names(nreg)[length(nreg)]
TOT_GAP <- strrep(" ", 4)
tot <- data.frame(label = rows$label, band = rows$band, half = rows$half,
                  region = factor(last_reg, levels = lev),
                  x = max(cols$col[cols$region == last_reg]),
                  # trimws because format() pads every number to the widest one's width,
                  # and those leading spaces make left-justified text read as right-aligned.
                  # trimws first: format() pads every number to the widest one's width, and
                  # those leading spaces make left-justified text read as right-aligned.
                  text = paste0(TOT_GAP,
                                trimws(format(as.integer(rows$n_runs), big.mark = ","))))

# Colour is carried as a literal hex per cell, so the scale is the identity: no palette is
# chosen here and none can drift from 02_reference's.
fills <- sort(unique(d$fill))

p <- ggplot(d, aes(x = col, y = label, fill = fill)) +
  geom_tile(width = 1, height = 0.82) +
  scale_fill_manual(values = setNames(fills, fills), guide = "none") +
  # space = "free" is the whole reason this is R: a region's panel is as wide as the runs it
  # holds, and the gaps between panels are the block separation that was hand-placed before.
  # ggh4x's nested strips give the Pathogens / Residents bracket for free, as an outer
  # strip over the band strips. Hand-drawn in the matplotlib version, this was a bracket
  # placed in absolute inches that had to be recomputed whenever a band changed height.
  # NO TINTED BOX behind the labels. It was carrying the band's modal matrix hue, and now
  # that every drawn row is a wheat pathogen every box was the same colour, so it encoded
  # nothing and competed with the cells for attention.
  geom_blank(data = pins, aes(x = x, y = label), inherit.aes = FALSE) +
  # The grid's border, drawn per band x region rather than by the theme, so the label
  # column stays unboxed. -Inf/Inf makes each rect its own panel exactly.
  geom_rect(data = frames, inherit.aes = FALSE, fill = NA, colour = "#8f8f8f",
            linewidth = 0.45, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_text(data = tot, aes(x = x, y = label, label = text), inherit.aes = FALSE,
            hjust = 0, size = 3.6, colour = "#3a3a3a") +
  # The category sits at the VERTICAL MIDDLE of its box, beside the taxa rather than above
  # them, and the taxa are LEFT-justified in a column of their own.
  geom_text(data = bands_df, aes(label = text, y = ymid), inherit.aes = FALSE,
            x = LAB_X_CAT, hjust = 0, vjust = 0.5, size = 3.6,
            fontface = "bold", colour = "#1a1a1a") +
  geom_text(data = labs, aes(x = x, y = label, label = text, fontface = face),
            inherit.aes = FALSE, hjust = 0, size = 3.6, colour = "#1a1a1a") +
  scale_discrete_identity(aesthetics = "fontface") +
  # THE Pathogens / Residents BRACKET IS GONE. Residents lost their rows, so every taxon
  # row left is a pathogen and the outer strip was a box enclosing the whole panel to say
  # so. One nesting layer remains, the band, and it is empty by the labeller because the
  # band name lives in the label box beside the taxa.
  facet_nested(band ~ region, scales = "free", space = "free_y", switch = "y",
               nest_line = element_blank(),
               strip = strip_nested(size = "variable", by_layer_y = TRUE,
                                    background_y = elem_list_rect(fill = NA, colour = NA,
                                                                  linewidth = 0)),
               labeller = labeller(band = function(x) rep("", length(x)))) +
  force_panelsizes(cols = widths) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  coord_cartesian(clip = "off") +
  # No axis title. The region strips already say what the x axis is, and the run count
  # belongs in the caption where it is computed rather than restated under the panel.
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    panel.spacing.x  = unit(5, "pt"),
    panel.spacing.y  = unit(7, "pt"),
    # NO panel.border. It is a theme setting and therefore applies to every panel, which
    # boxed the label column as if it were a region of data. The grid panels get their
    # border from `frames` below, drawn only where there are cells.
    panel.border     = element_blank(),
    axis.text.y      = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks       = element_blank(),
    strip.text.x     = element_text(size = 11, face = "bold", colour = "#1a1a1a", angle = 0,
                                    margin = margin(b = 4)),
    # The band strip is emptied by the LABELLER, not by element_blank(): ggh4x sizes the
    # nested strips from the element, and blanking it leaves the size NULL and the build
    # fails. The band name is in the label box now; only the Pathogens/Residents strip shows.
    strip.text.y.left = element_text(size = 10, face = "bold", colour = "#3a3a3a", angle = 90,
                                     margin = margin(r = 3, l = 3)),
    strip.placement  = "outside",
    plot.margin      = margin(6, 54, 6, 6))

# Width is set from the content: the widest region strip label plus the row labels plus the
# grid. 10.2 in clipped Oceania's panel off the right edge, which is the kind of thing a
# fixed size does silently.
H <- 0.75 + 0.17 * nrow(rows) + 0.20 * nlevels(rows$band)
ragg::agg_png(file.path(HERE, "samples.png"), width = 11.4, height = H,
              units = "in", res = 600, background = "white")
print(p)
invisible(dev.off())

# ── legend ────────────────────────────────────────────────────────────────────
# THREE keys, because the panel uses colour for three different things and a legend that
# showed only one of them was telling the reader a third of the story:
#
#   hue       which category of the reference matrix the taxon is
#   treatment whether the study that generated the run declared it
#   ramp      the per-run totals, one scale per half, NOT comparable with each other
#
# Every colour and every name is read from the tables, so a category that stops occurring
# stops appearing and none of this can drift from what was drawn.
ramps <- rd("samples_ramps.tsv")
ramps$k <- as.integer(ramps$k); ramps$max <- as.integer(ramps$max)

# The matrix-category keys, PLUS the study-target marker, which is a different kind of
# fact and so is named as what the study did rather than as what an organism is.
cats <- unique(rows[!rows$taxon %in% "__total_residents__", c("category", "colour")])
cats <- cats[order(match(cats$colour, c(unique(rows$colour)))), ]
kd1 <- data.frame(k = factor(cats$category, levels = cats$category), x = 1, y = 1)
g1 <- ggplot(kd1, aes(x, y, fill = k)) + geom_tile() +
  scale_fill_manual(values = setNames(cats$colour, cats$category), name = NULL,
                    guide = guide_legend(nrow = 1, keywidth = unit(11, "pt"),
                                         keyheight = unit(11, "pt")))

wp  <- rows[rows$half == "Pathogens", ][1, ]
trt <- c("declared by the study"   = unique(cells$fill[cells$state == "declared"])[1],
         "not declared"            = wp$colour,
         "declaration not curated" = wp$colour_unknown)
trt <- trt[c("declared", "undeclared", "unknown") %in% unique(cells$state)]
kd2 <- data.frame(k = factor(names(trt), levels = names(trt)), x = 1, y = 1)
g2 <- ggplot(kd2, aes(x, y, fill = k)) + geom_tile() +
  scale_fill_manual(values = trt, name = NULL,
                    guide = guide_legend(nrow = 1, keywidth = unit(11, "pt"),
                                         keyheight = unit(11, "pt")))

# The two ramps are drawn as the discrete steps they are: the counts are small integers, so
# a continuous bar would imply a precision the data does not have.
# The range goes in the LABEL, not on the tiles. Drawn on the tiles, the end numbers were
# offset vertically to clear the bar and landed in the other half's row.
rng <- setNames(sprintf("%s per run (0\u2013%d)", tolower(unique(ramps$half)),
                        tapply(ramps$max, ramps$half, max)[unique(ramps$half)]),
                unique(ramps$half))
g3 <- ggplot(ramps, aes(x = k, y = half, fill = fill)) +
  geom_tile(colour = "#d8d8d8", linewidth = 0.2) +
  scale_fill_manual(values = setNames(ramps$fill, ramps$fill), guide = "none") +
  scale_y_discrete(labels = function(x) rng[x]) +
  theme_void(base_size = 11) +
  theme(axis.text.y = element_text(size = 9.5, colour = "#1a1a1a",
                                   hjust = 1, margin = margin(r = 5)),
        plot.margin = margin(2, 2, 2, 2))

# ONE ROW, Leon's call. Every guide is nrow = 1 above and the three blocks sit side by side
# here, so the legend is a single band under the panel rather than a stacked block beside
# it. That costs width and buys height, which is the scarce axis: the composite puts four
# detection-space panels above this one on a single A4 page.
pick <- function(p) cowplot::get_plot_component(
  p + theme_void() + theme(legend.position = "bottom",
                           # the legend is an annotation, not a heading: regular weight
                           legend.text = element_text(size = 9.5, colour = "#1a1a1a",
                                                      margin = margin(l = 4, r = 26)),
                           legend.key.spacing.x = unit(2, "pt"),
                           # a zero left margin clipped the first key's swatch off the
                           # edge of its grid cell, leaving the first label unexplained
                           legend.margin = margin(0, 0, 0, 6)),
  "guide-box-bottom", return_all = TRUE)

# Widths are proportional to the TEXT each block holds, so the three sit next to each other
# rather than floating in the middle of equal thirds.
lg <- cowplot::plot_grid(pick(g2), pick(g1), g3, nrow = 1, rel_widths = c(1.15, 1.5, 1.2))
ragg::agg_png(file.path(HERE, "samples_legend.png"), width = 11.4, height = 0.32,
              units = "in", res = 600, background = "white")
print(lg)
invisible(dev.off())

cat(sprintf("wrote %s/samples.png\n", HERE))
cat(sprintf("  %d cells, %d rows over %d bands, %d region panels, %d runs\n",
            nrow(d), nrow(rows), nlevels(rows$band), nlevels(cols$region), nrow(cols)))
