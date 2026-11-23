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
LAB_X_CAT  <- 0.02                     # category name, left edge
LAB_X_TAXA <- 0.34                     # taxa column, left-justified
labs <- data.frame(x = LAB_X_TAXA, label = rows$label, band = rows$band, half = rows$half,
                   region = factor(LABEL_PANEL, levels = lev),
                   text = as.character(rows$label),
                   face = ifelse(rows$italic == "1", "italic", "plain"))
pins <- data.frame(x = c(0, 1), label = rows$label[1], band = rows$band[1],
                   half = rows$half[1], region = factor(LABEL_PANEL, levels = lev))
# The BAND NAME lives in the label box too, not in a strip beside it, so one box carries the
# category and the taxa it contains. One row per band, placed at the box's left edge and
# centred on the band's rows; the box is tinted with the band's prevailing matrix colour.
# The category is placed on the band's MIDDLE row, because a discrete y scale has no
# "centre" to address: each band panel holds only its own labels, so the midpoint is the
# label at ceiling(k/2). With an even row count it lands a half-row high, which is close
# enough and is stable.
# Deduplicate on HALF AND BAND, not band alone: "Total" occurs once per half, and keying
# on the name alone gave the second one no box and no label.
bl <- rows[!duplicated(rows[, c("half", "band")]), ]
mid <- mapply(function(h, b) {
  lb <- rows$label[rows$half == h & rows$band == b]
  as.character(lb[ceiling(length(lb) / 2)])
}, bl$half, bl$band)
bands_df <- data.frame(band = factor(bl$band, levels = levels(rows$band)),
                       half = factor(bl$half, levels = levels(rows$half)),
                       region = factor(LABEL_PANEL, levels = lev),
                       fill = bl$band_fill, text = as.character(bl$band),
                       ymid = factor(mid, levels = levels(rows$label)))
LAB_SHARE <- 0.21                      # of the drawn width, enough for the longest binomial
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
fills <- sort(unique(c(d$fill, rows$band_fill)))

p <- ggplot(d, aes(x = col, y = label, fill = fill)) +
  geom_tile(width = 1, height = 0.82) +
  scale_fill_manual(values = setNames(fills, fills), guide = "none") +
  # space = "free" is the whole reason this is R: a region's panel is as wide as the runs it
  # holds, and the gaps between panels are the block separation that was hand-placed before.
  # ggh4x's nested strips give the Pathogens / Residents bracket for free, as an outer
  # strip over the band strips. Hand-drawn in the matplotlib version, this was a bracket
  # placed in absolute inches that had to be recomputed whenever a band changed height.
  geom_rect(data = bands_df, aes(fill = fill), inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_blank(data = pins, aes(x = x, y = label), inherit.aes = FALSE) +
  geom_text(data = tot, aes(x = x, y = label, label = text), inherit.aes = FALSE,
            hjust = 0, size = 2.6, colour = "#4a4a4a") +
  # The category sits at the VERTICAL MIDDLE of its box, beside the taxa rather than above
  # them, and the taxa are LEFT-justified in a column of their own.
  geom_text(data = bands_df, aes(label = text, y = ymid), inherit.aes = FALSE,
            x = LAB_X_CAT, hjust = 0, vjust = 0.5, size = 2.7, colour = "#4a4a4a") +
  geom_text(data = labs, aes(x = x, y = label, label = text, fontface = face),
            inherit.aes = FALSE, hjust = 0, size = 2.75, colour = "#1a1a1a") +
  scale_discrete_identity(aesthetics = "fontface") +
  facet_nested(half + band ~ region, scales = "free", space = "free_y", switch = "y",
               nest_line = element_blank(),
               # The OUTER strip gets a box, so Pathogens and Residents each enclose the
               # band rows they cover instead of being a word floating beside them. The
               # inner strip is empty (the band name is in the label box) and is left
               # unpainted. elem_list_rect applies per nesting layer, innermost first.
               strip = strip_nested(
                 size = "variable",
                 # by_layer_y is required: without it elem_list_rect recycles across
                 # every strip in draw order, so the two colours alternate down the panel
                 # instead of distinguishing the outer layer from the inner one. Verified
                 # by rendering a two-colour test rather than reasoning about it.
                 by_layer_y = TRUE,
                 # Element 1 is the OUTER layer (leftmost under switch="y"), which is the
                 # one carrying the Pathogens / Residents text; element 2 is the empty band
                 # strip. Painted the other way round the box landed beside the words
                 # instead of round them.
                 background_y = elem_list_rect(fill = c("#ffffff", NA),
                                               colour = c("#8f8f8f", NA),
                                               linewidth = c(0.45, 0))),
               labeller = labeller(band = function(x) rep("", length(x)),
                                   half = identity)) +
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
    panel.border     = element_rect(fill = NA, colour = "#8f8f8f", linewidth = 0.45),
    axis.text.y      = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks       = element_blank(),
    strip.text.x     = element_text(size = 8, colour = "#1a1a1a", angle = 0,
                                    margin = margin(b = 4)),
    # The band strip is emptied by the LABELLER, not by element_blank(): ggh4x sizes the
    # nested strips from the element, and blanking it leaves the size NULL and the build
    # fails. The band name is in the label box now; only the Pathogens/Residents strip shows.
    strip.text.y.left = element_text(size = 7.5, colour = "#4a4a4a", angle = 90,
                                     margin = margin(r = 3, l = 3)),
    strip.placement  = "outside",
    plot.margin      = margin(6, 54, 6, 6))

# Width is set from the content: the widest region strip label plus the row labels plus the
# grid. 10.2 in clipped Oceania's panel off the right edge, which is the kind of thing a
# fixed size does silently.
H <- 1.05 + 0.20 * nrow(rows) + 0.22 * nlevels(rows$band)
ragg::agg_png(file.path(HERE, "samples.png"), width = 11.4, height = H,
              units = "in", res = 600, background = "white")
print(p)
invisible(dev.off())

# ── legend ────────────────────────────────────────────────────────────────────
# Built as a REAL ggplot guide and extracted with cowplot, not as hand-placed text at
# invented coordinates. The hand-placed version drifted as soon as a label changed length,
# because nothing was measuring anything; a guide lays itself out from the strings it is
# given.
#
# Two keys and the split is the point: the TREATMENT of a colour says who knew, the HUE says
# which reference-matrix category the taxon is. Both values come from the tables.
wp  <- rows[rows$half == "Pathogens", ][1, ]
# Only the states PRESENT in the cells get a key. Once the declarations were curated the
# "not curated" state stopped occurring, and a key for a colour the panel never shows is a
# reader looking for something that is not there.
key <- c("declared by the study"   = unique(cells$fill[cells$state == "declared"])[1],
         "not declared"            = wp$colour,
         "declaration not curated" = wp$colour_unknown)
present <- c("declared", "undeclared", "unknown") %in% unique(cells$state)
key <- key[present]
kd  <- data.frame(k = factor(names(key), levels = names(key)), x = 1, y = 1)

g <- cowplot::get_plot_component(
  ggplot(kd, aes(x, y, fill = k)) + geom_tile() +
    scale_fill_manual(values = key, name = NULL,
                      guide = guide_legend(nrow = 1, byrow = TRUE,
                                           keywidth = unit(11, "pt"),
                                           keyheight = unit(8, "pt"))) +
    theme_void(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7.5, colour = "#1a1a1a",
                                     margin = margin(l = 3, r = 10)),
          legend.key.spacing.x = unit(2, "pt")),
  "guide-box-bottom", return_all = TRUE)

ragg::agg_png(file.path(HERE, "samples_legend.png"), width = 4.6, height = 0.26,
              units = "in", res = 600, background = "white")
grid.newpage(); grid.draw(g)
invisible(dev.off())

cat(sprintf("wrote %s/samples.png\n", HERE))
cat(sprintf("  %d cells, %d rows over %d bands, %d region panels, %d runs\n",
            nrow(d), nrow(rows), nlevels(rows$band), nlevels(cols$region), nrow(cols)))
