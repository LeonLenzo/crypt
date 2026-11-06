#!/usr/bin/env Rscript
# map_cohort.R — where the wheat field cohort was collected.
#
#   Rscript release/01_cohort/scripts/map_cohort.R
# Inputs : figures/cohort_map_points.tsv, cohort_map_summary.tsv  (prep_cohort_map.py)
# Output : figures/cohort_map.{png,pdf}
#          figures/cohort_map_legend.{png,pdf}
#
# THE LEGEND IS A SEPARATE FILE, as it is for the network figure this sits beside. The
# two are laid out together, and a key baked into one panel cannot be placed to suit
# both. LEGEND_ON_MAP puts it back for a quick look.
#
# This replaces map_crop_origins.R, which drew four cereals. That figure needed colour
# for crop, and needed bar panels beside it, because a map of four crops with wildly
# different sampling shapes reads as "mostly wheat" on its own: rice was 2,784 samples
# at 6 localities, wheat 1,251 at 378. None of that applies now. The cohort is wheat,
# so colour carries nothing, the comparison panels have nothing to compare, and what is
# left is the thing the panels existed to let you see anyway — that this cohort is broad
# rather than deep, 1,098 runs spread over 441 places.
#
# POINT SIZE IS LOG, not area, carried over from the crop map and for the same reason.
# Counts here run 1 to 30, which is tame next to the 1 to 1,396 that forced the decision
# originally, but the argument holds: most places contributed one or two runs and they
# are the substance of the figure, so they have to stay visible. The cost belongs in any
# caption: a circle twice the area is NOT twice the samples. It is a magnitude cue.
#
# TWO COLOURS, ONE SHAPE. Green is a named place, blue is a country centroid, for the
# 289 runs whose papers record a country and nothing finer. Diamonds did this job first
# and were swapped out: at these sizes a diamond and a circle are hard to tell apart in
# a cluster, and the difference matters most exactly where the map is densest. Both
# colours are the network figure's, so panels A and B speak the same language.
#
# THE BIGGEST MARKS ARE THE LEAST PRECISE, which is the opposite of what a reader
# expects, so it has to be said rather than designed around. South Africa's 65-run
# diamond is the largest symbol on the map and it means "somewhere in South Africa".
# Both series share one size scale so the quantities stay comparable.
#
# A DIAMOND IS NOT A FIELD and the caption has to say so. Those runs were chased hard —
# three papers' supplements read end to end — and the resolution is the literature's,
# not ours. Drawing them as circles would claim a precision that does not exist;
# leaving them off would hide a fifth of the cohort.
#
# WHAT IS STILL MISSING. 32 runs name no place at all. That is in the summary table.

suppressMessages({
  library(ggplot2); library(maps); library(ragg); library(cowplot)
})

here    <- "release/01_cohort"
pts     <- read.delim(file.path(here, "figures/cohort_map_points.tsv"),
                      stringsAsFactors = FALSE)
summary <- read.delim(file.path(here, "figures/cohort_map_summary.tsv"),
                      stringsAsFactors = FALSE)
get <- function(name) summary$value[summary$measure == name]

LEGEND_ON_MAP <- FALSE

LAND   <- "grey92"      # filled, not outlined: the points are the subject
BORDER <- "white"
POINT   <- "#4a9b34"    # the field-wheat green the network figure uses
COUNTRY <- "#2a78d6"    # the BioProject blue from the same figure
OUTLINE <- "black"
# Where the magnifier sits. Its proportions are not free: annotation_custom stretches
# the grob to the box, and the inset holds its own coord_fixed, so a box of the wrong
# shape letterboxes the content and leaves white bands inside the frame. At this crop
# the content wants a 1.41 display ratio, which with the main map's own 1.25 fixes the
# box at about 104 degrees of longitude per 59 of latitude. Reaching further east over
# Russia was tried and is what produced the bands. Sliding the box east keeps the shape
# and so keeps it filled; it was moved until it sat clear of the Americas without
# reaching the Pakistani or Ethiopian sites on the other side.
INSET_BOX <- c(-58, 46, 19, 78)

world <- map_data("world")
world <- world[world$region != "Antarctica", ]

p <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = LAND, colour = BORDER, linewidth = 0.15) +
  geom_point(data = pts, aes(lon, lat, size = runs, fill = placement),
             shape = 21, colour = OUTLINE, stroke = 0.4, alpha = 0.78) +
  scale_fill_manual(name = NULL, values = c("Locality" = POINT, "Country" = COUNTRY),
                    labels = c("Locality" = "a named place",
                               "Country" = "country only")) +
  # Breaks reach 65 because the largest marks on the map are country diamonds, not
  # localities: diamonds carry 1-65 runs and circles only 1-30. A key that stopped at
  # 30 left the biggest symbols on the figure unexplained.
  scale_size_continuous(name = "runs", trans = "log10", limits = c(1, 65),
                        range = c(1.2, 7), breaks = c(1, 5, 20, 65)) +
  coord_fixed(1.25, xlim = c(-170, 180), ylim = c(-56, 78), expand = FALSE) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    legend.position  = if (LEGEND_ON_MAP) c(0.08, 0.34) else "none",
    legend.spacing.y = unit(2, "pt")) +
  guides(shape = guide_legend(override.aes = list(size = 4, fill = POINT), order = 1),
         size  = guide_legend(override.aes = list(shape = 21, fill = POINT), order = 2))

# INSET: EUROPE, NORTH AFRICA AND THE NEAR EAST. 289 of the 460 points and 852 of the
# runs sit in this window, so at world scale most of the cohort is one blot.
#
# The window is cropped to the data, not to the region: south to Tenerife at 28.5N and
# east to Turkey at 35.5E, which are the real corners of what is in it. Rounding out to
# a tidy "Europe and the Mediterranean" box spent a third of the panel on empty Sahara
# and Atlantic and shrank everything worth seeing.
#
# IT SITS OVER ITSELF, extended west across the Atlantic. That placement only works
# because the window was widened until it CONTAINS everything the placement box covers:
# Morocco, the Canaries and Israel were each hidden by earlier, narrower versions, and
# a magnifier that deletes three real sites to show the rest is worse than none. The
# box was then checked against every point on the map, and hides nothing it does not
# itself redraw, larger. It was grown to the largest size that keeps that true, east
# across Russia to 68E and west to 72W: six degrees further east and two Pakistani
# sites vanish behind it, one degree further west and three North American ones do.
#
# An earlier version floated this over the Indian Ocean, picked by searching for empty
# canvas. It read as a second, unexplained map. Over its own ground it reads as a
# magnifier.
EU <- c(-19, 38, 27, 58)
inset <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = LAND, colour = BORDER, linewidth = 0.15) +
  geom_point(data = pts, aes(lon, lat, size = runs, fill = placement),
             shape = 21, colour = OUTLINE, stroke = 0.4, alpha = 0.78) +
  scale_fill_manual(values = c("Locality" = POINT, "Country" = COUNTRY),
                    guide = "none") +
  scale_size_continuous(trans = "log10", limits = c(1, 65), range = c(1.2, 7),
                        guide = "none") +
  coord_fixed(1.3, xlim = EU[1:2], ylim = EU[3:4], expand = FALSE) +
  theme_void() +
  # The frame has to be on plot.background, not panel.background. On the panel it is
  # drawn inside a transparent margin, so the world map behind shows through around it
  # and the coastlines appear to cross the inset. Zero margin puts the panel edge and
  # the frame in the same place.
  theme(plot.background  = element_rect(fill = "white", colour = "black",
                                        linewidth = 1.6),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin      = margin(0, 0, 0, 0))

# Drawn after the world layer but with the points still visible around it; the box is
# opaque, so whatever it covers is replaced by the magnified version inside it.
p <- p +
  annotation_custom(ggplotGrob(inset), xmin = INSET_BOX[1], xmax = INSET_BOX[2],
                    ymin = INSET_BOX[3], ymax = INSET_BOX[4]) +
  # The frame is drawn here, not inside the inset. annotation_custom clips the grob to
  # its box, so a border on the inset's own plot.background loses the outer half of its
  # stroke and renders as a hairline whatever width it is set to.
  annotate("rect", xmin = INSET_BOX[1], xmax = INSET_BOX[2],
           ymin = INSET_BOX[3], ymax = INSET_BOX[4],
           fill = NA, colour = "black", linewidth = 0.9)

agg_png(file.path(here, "figures/cohort_map.png"), width = 11, height = 6,
        units = "in", res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(here, "figures/cohort_map.pdf"), width = 11, height = 6)
print(p); invisible(dev.off())

# The key on its own. Built by drawing the same scales over an empty panel and keeping
# only the guide box, so it cannot drift from what the map actually drew.
key <- ggplot(pts, aes(lon, lat, size = runs, fill = placement)) +
  geom_point(shape = 21, colour = OUTLINE, stroke = 0.4, alpha = 0.78) +
  scale_fill_manual(name = NULL, values = c("Locality" = POINT, "Country" = COUNTRY),
                    labels = c("Locality" = "Locality",
                               "Country" = "Country")) +
  scale_size_continuous(name = "runs", trans = "log10", limits = c(1, 65),
                        range = c(1.2, 7), breaks = c(1, 5, 20, 65)) +
  guides(shape = guide_legend(override.aes = list(size = 4, fill = POINT), order = 1),
         size  = guide_legend(override.aes = list(shape = 21, fill = POINT), order = 2)) +
  theme_minimal(base_size = 16) +
  theme(legend.background = element_rect(fill = "white", colour = NA))
# get_plot_component can return a list of candidate grobs; take the first real one.
legend <- cowplot::get_plot_component(key, "guide-box-right", return_all = TRUE)
if (inherits(legend, "list")) legend <- legend[[which(
  vapply(legend, function(g) inherits(g, "gtable") || inherits(g, "grob"), logical(1)))[1]]]

agg_png(file.path(here, "figures/cohort_map_legend.png"), width = 2.4, height = 3.2,
        units = "in", res = 300, background = "white")
grid::grid.newpage(); grid::grid.draw(legend); invisible(dev.off())
cairo_pdf(file.path(here, "figures/cohort_map_legend.pdf"), width = 2.4, height = 3.2)
grid::grid.newpage(); grid::grid.draw(legend); invisible(dev.off())

cat(sprintf("-> cohort_map.{png,pdf}  %s points; %s runs at a locality, %s at a country centroid, of %s field-wheat runs\n",
            nrow(pts), get("drawn at a locality"),
            get("drawn at a country centroid"), get("field-wheat runs in cohort")))
cat(sprintf("   %s runs name no place at all and are NOT drawn\n",
            get("not drawn: no usable place")))
