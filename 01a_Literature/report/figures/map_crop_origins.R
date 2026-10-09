#!/usr/bin/env Rscript
# map_crop_origins.R — where the field cereal cohort was collected, by crop.
#
# The 02_literature map (map_sample_origins.R) put one circle per COUNTRY and filled it by
# collection period, which works because period is a single ordered value per country. Crop is
# categorical and several crops share a locality, so this plots one point per locality x crop
# at geocoded coordinates instead. Country circles cannot show crop composition without pie
# charts, and pies at this size are unreadable.
#
# The structure this figure has to survive: the crops are not comparable in shape. Rice is
# 2,784 samples across 6 localities and wheat is 1,251 across 378, because the rice data are
# large GWAS panels grown at one farm while the wheat data are rust surveys that sample
# broadly. A map alone therefore reads as "mostly wheat" and a bar chart alone reads as
# "mostly rice". Both panels are present for that reason, and the localities panel exists
# because it is the only one that shows the difference.
#
# Point size is on a LOG scale, not an area scale (leon, 2026-10-09: "much less severe,
# larger at the bottom smaller at the top"). Sample counts span 1 to 1,396, so area scaling
# gives a 37-fold radius ratio and three big circles swamp the 378 single-digit wheat
# localities that are the point of the panel. log10 across a radius range of 1.8 to 7 puts
# that ratio at under 4: a one-sample locality is still a visible dot and a 1,396-sample one
# is merely the biggest.
#
# The cost, and it belongs in any caption: a log size scale is NOT area-proportional, so a
# circle twice the area does not mean twice the samples. It is a magnitude cue, not a
# quantity. The bar panels carry the actual numbers.
#
# The Europe inset is gone (leon, 2026-10-09): once the points were log-scaled, borderless
# and at alpha 0.6, the European cluster resolves at world scale and the inset was repeating
# the main panel rather than adding to it.
#
# Crops under 100 samples are dropped by prep_crop_map.py, not here, so the bar panels and the
# map always agree about what exists. What was dropped is in crop_totals.tsv with an
# `included` flag; as of 2026-10-08 that is barley (3) and the small cereals (53).
#
#   Rscript 01a_Literature/report/figures/map_crop_origins.R
# Inputs : data/gold/crop_map.tsv, data/gold/crop_totals.tsv  (prep_crop_map.py)
# Output : report/figures/map_crop_origins.{png,pdf}
#
# Panels, left to right under the map: samples per crop with a REVERSED x axis, localities
# per crop, then collection year as a stacked area. The first two are back to back so the
# crop labels sit once, between them, and the two quantities read as one comparison instead
# of two charts that happen to share categories.

# --panels writes each panel to its own file as well as the composite, so the figure can be
# laid out by hand instead of by patchwork.
PANELS <- "--panels" %in% commandArgs(TRUE)

suppressMessages({
  library(ggplot2); library(maps); library(dplyr); library(scales)
  library(patchwork); library(ragg)
})

here <- "01a_Literature"
pts  <- read.delim(file.path(here, "data/gold/crop_map.tsv"), stringsAsFactors = FALSE)
tot  <- read.delim(file.path(here, "data/gold/crop_totals.tsv"), stringsAsFactors = FALSE)
yr   <- read.delim(file.path(here, "data/gold/crop_year.tsv"), stringsAsFactors = FALSE)

tot <- tot[tot$included == "yes", ]
ORDER <- c("Wheat", "Maize", "Rice", "Sorghum", "Barley", "Other cereal",
           "Cereal, host unresolved")
ORDER <- ORDER[ORDER %in% unique(c(pts$crop, tot$crop))]
pts$crop <- factor(pts$crop, levels = ORDER)
tot$crop <- factor(tot$crop, levels = ORDER)
yr$crop  <- factor(yr$crop,  levels = ORDER)

# Muted, semantic palette (leon, 2026-10-09: keep the hues, lose the primary-school
# saturation). The hue assignments are the crop's own: wheat grown cold and blue, maize
# kernels gold, rice harvested green, sorghum seed a rust orange.
#
# Chosen by measurement, not eye. Against the Okabe-Ito set it replaces, mean chroma drops
# from 61 to 43 - which is the "polished" part - while the worst-separated PAIR under
# simulated deuteranopia improves from dE 11 to 28. Okabe-Ito is colourblind-safe as an
# eight-colour set, but the four entries this figure used happened to put maize (#E69F00)
# and sorghum (#D55E00) in nearly the same place once red-green vision is removed.
fills <- c("Wheat"                   = "#3A6EA5",   # steel blue
           "Maize"                   = "#D9A441",   # old gold
           "Rice"                    = "#45896A",   # sage green
           "Sorghum"                 = "#A9512C",   # burnt sienna
           "Barley"                  = "#8C6D8F",
           "Other cereal"            = "#8D9EA9",
           "Cereal, host unresolved" = "grey62")

OCEAN <- "white"      # also the internal border colour; see map_layers

# Chart panels carry drawn axes rather than gridlines (leon, 2026-10-09). theme_minimal has
# no axis lines at all, so they have to be added back explicitly; the grid is removed in the
# same block so the two can never both be on.
# Lines, ticks and text are all black so the axis reads as one object (leon, 2026-10-09).
# theme_minimal sets axis text to grey30, so the text is set here too rather than leaving
# black rules against grey numerals.
# Two breaks, zero and the next round number down from the maximum. pretty() still returned
# four in a square panel, which collided into "01,00020003000" in the composite even though
# it read fine in the 5x5 standalone panel.
two_breaks <- function(x) {
  m <- max(x, na.rm = TRUE)
  u <- 10 ^ floor(log10(m))
  c(0, floor(m / u) * u)
}

axes <- theme(panel.grid = element_blank(),
              axis.line = element_line(colour = "black", linewidth = 0.6),
              axis.ticks = element_line(colour = "black", linewidth = 0.6),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              axis.ticks.length = unit(4, "pt"),
              # The last x-axis label sits at the panel edge and is half a label wide, so
              # without a right margin "3,000" renders as "3,00".
              plot.margin = margin(6, 16, 6, 6))

base <- theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        plot.title       = element_text(size = 15, margin = margin(b = 8)))

world <- map_data("world")
world <- world[world$region != "Antarctica", ]

# Rows with precision == "country" carry no coordinates: they name a country and nothing
# finer, so prep_crop_map.py leaves the placement to the figure. Join the centroid of each
# country's LARGEST polygon, not the midpoint of its bounding box, which would put the USA
# in the Pacific and France in the Atlantic.
cent <- world %>%
  group_by(region, group) %>%
  summarise(long = mean(long), lat = mean(lat), npts = n(), .groups = "drop") %>%
  group_by(region) %>% slice_max(npts, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(region, clon = long, clat = lat)

# "Zimbabwae" is not a variant spelling, it is a typo in the submitted BioSample metadata,
# and it is mapped here rather than corrected upstream because the archive string is what
# the curation records. Three samples.
TO_MAPS <- c("United Kingdom" = "UK", "USA" = "USA", "Czechia" = "Czech Republic",
             "Zimbabwae" = "Zimbabwe", "The Netherlands" = "Netherlands")
pts$region <- ifelse(pts$country %in% names(TO_MAPS), TO_MAPS[pts$country], pts$country)
pts <- pts %>% left_join(cent, by = "region") %>%
  mutate(lat = ifelse(precision == "country", clat, lat),
         lon = ifelse(precision == "country", clon, lon),
         placement = ifelse(precision == "country", "Country", "Locality"))

unplaced_ctry <- pts %>% filter(precision == "country", is.na(lat))
if (nrow(unplaced_ctry)) message("country not on the basemap: ",
                                 paste(unique(unplaced_ctry$country), collapse = ", "))
pts <- pts %>% filter(!is.na(lat))
pts$placement <- factor(pts$placement, levels = c("Locality", "Country"))

pts <- pts %>% arrange(desc(n))

n_cer <- sum(tot$runs); n_loc <- length(unique(pts$location)); n_pl <- sum(pts$n)

SIZE_RANGE  <- c(1.8, 7)
SIZE_BREAKS <- c(1, 10, 100, 1000)

# Country outlines are drawn in the OCEAN colour rather than a grey (leon, 2026-10-09).
# Borders then read as hairline separations in the landmass instead of as drawn lines, and
# nothing in the basemap competes with the points.
map_layers <- function(xlim, ylim, ratio, srange) {
  list(
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "grey92", colour = OCEAN, linewidth = 0.35),
    # Shape, not colour, separates a real locality from a country centroid: a centroid is a
    # statement about which COUNTRY, drawn at a place nobody sampled, and a reader must not
    # mistake the diamond in the middle of France for a field site.
    geom_point(data = pts, aes(lon, lat, fill = crop, size = n, shape = placement),
               stroke = 0, alpha = 0.6),
    scale_shape_manual(values = c("Locality" = 21, "Country" = 23),
                       name = "Resolution", drop = FALSE),
    scale_size(transform = "log10", range = srange, breaks = SIZE_BREAKS,
               labels = comma, name = "Samples"),
    scale_fill_manual(values = fills, name = "Crop", drop = FALSE),
    coord_fixed(ratio, xlim = xlim, ylim = ylim, expand = FALSE)
  )
}

p_map <- ggplot() +
  # Full extent, not a crop to where the points are (leon, 2026-10-09). The old
  # c(-165, 180) x c(-48, 72) window cut the top off Greenland and the bottom off southern
  # Chile and New Zealand, which reads as a rendering fault. Whitespace above and below is
  # the acceptable cost. Antarctica is dropped from the basemap, not clipped by the window.
  map_layers(c(-180, 180), c(-57, 84), 1.32, SIZE_RANGE) +
  guides(fill  = guide_legend(override.aes = list(size = 5, shape = 21), order = 1, ncol = 1),
         shape = guide_legend(override.aes = list(size = 5, fill = "grey35"), order = 2),
         size  = guide_legend(override.aes = list(fill = "grey35", shape = 21), order = 3)) +
  base +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(),
        # Inside the map, over the empty South Pacific: the right margin now belongs to the
        # chart column, and a legend there would push the map into a strip.
        legend.position = "inside", legend.position.inside = c(0.012, 0.30),
        legend.justification = c(0, 0.5), legend.box = "vertical",
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.spacing.y = unit(2, "pt"))

# No bar outlines (leon, 2026-10-09): the fills are muted enough to hold their own shape and
# the black keylines were the last of the primary-school look.
#
# All three panels are horizontal bars reading left to right. The samples axis was reversed
# while b and c sat side by side back-to-back; stacked down the right-hand column there is
# nothing to mirror against, so a reversed axis would just be a second direction for a reader
# to track.
bar <- function(df, xvar, xlab) {
  ggplot(df, aes(.data[[xvar]], crop, fill = crop)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = fills, guide = "none") +
    # No x axis on b and c. Every bar carries its exact count, so the axis was repeating
    # the labels; in a square panel the two together collided into "01,00020003000" and
    # truncated the values to "2,". The axis TITLE stays, because it is what names the
    # quantity.
    # 0.34 was not enough headroom: the label sits outside the bar, and in a square panel
    # "2,784" is wider than the space left past the longest bar, so it clipped to "2,".
    scale_x_continuous(expand = expansion(c(0, 0.08)), labels = comma) +
    scale_y_discrete(limits = rev(ORDER)) +
    labs(x = xlab, y = NULL) +
    base + axes + theme(aspect.ratio = 1)
}

p_runs <- bar(tot, "runs", "Samples")
p_locs <- bar(tot, "localities", "Distinct localities")

# Stacked BARS, not an area (leon, 2026-10-09). Collection year is a discrete count, and an
# area chart draws a slope between 2015 and 2016 that asserts samples were collected at every
# instant in between. Nothing was. Bars also let a year with no collection simply be absent
# instead of a V-shaped notch.
#
# Horizontal, like b and c, so all three panels in the column share an orientation and the
# reader changes axis once rather than three times. Earliest year at the top, so the column
# reads downward as time moves forward.
p_year <- ggplot(yr, aes(n, factor(year), fill = crop)) +
  geom_col(width = 0.76) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_x_continuous(labels = comma, expand = expansion(c(0, 0.04)),
                     breaks = function(x) pretty(x, 3)) +
  # Every second year is labelled; all 13 bars are still drawn. Thirteen labels in a square
  # panel collide, and the gap between pulses is legible from the bars themselves.
  scale_y_discrete(limits = rev(sort(unique(as.character(yr$year)))),
                   breaks = function(x) x[seq(length(x), 1, by = -2)]) +
  labs(x = "Samples", y = "Collection year") +
  base + axes + theme(axis.text.y = element_text(size = 12), aspect.ratio = 1)

# The right column is assembled FIRST and given its own heights. Written as
# `p_map | (a / b / c) + plot_layout(...)` the layout applied to the top-level row, which has
# one row, so the heights were silently ignored and the three charts collapsed into a strip.
right <- p_runs / p_locs / p_year + plot_layout(heights = c(1, 1, 2.6))

out <- (p_map | right) +
  plot_layout(widths = c(3.4, 1)) +
  # The title belongs to the FIGURE, not to panel (a): set on the map it collided with the
  # panel tag and rendered as "aDistribution of...".
  plot_annotation(title = "Distribution of Field Cereal RNA-seq Studies",
                  tag_levels = "a",
                  theme = theme(plot.title = element_text(
                    face = "bold", size = 23, margin = margin(b = 12)),
                    plot.background = element_rect(fill = "white", colour = NA))) &
  theme(plot.tag = element_text(face = "bold", size = 19),
        plot.tag.position = c(0, 1))

agg_png(file.path(here, "report/figures/map_crop_origins.png"), width = 17, height = 7.4,
        units = "in", res = 300, background = "white")
print(out); invisible(dev.off())
cairo_pdf(file.path(here, "report/figures/map_crop_origins.pdf"), width = 17, height = 7.4)
print(out); invisible(dev.off())
cat("wrote report/figures/map_crop_origins.{png,pdf}\n")

if (PANELS) {
  one <- function(g, stem, w, h) {
    agg_png(file.path(here, sprintf("report/figures/panels/%s.png", stem)), width = w,
            height = h, units = "in", res = 300, background = "white")
    print(g); invisible(dev.off())
    cairo_pdf(file.path(here, sprintf("report/figures/panels/%s.pdf", stem)),
              width = w, height = h)
    print(g); invisible(dev.off())
    cat(sprintf("  panels/%s.{png,pdf}  %gx%g in\n", stem, w, h))
  }
  dir.create(file.path(here, "report/figures/panels"), showWarnings = FALSE)
  one(p_map  + theme(legend.position = "right"), "a_map",        12, 6.2)
  one(p_runs, "b_samples",    5, 5)
  one(p_locs, "c_localities", 5, 5)
  one(p_year, "d_year",       5, 5)
}
cat(sprintf("  %d points, %s of %s samples placed, %d localities\n",
            nrow(pts), comma(n_pl), comma(n_cer), n_loc))
