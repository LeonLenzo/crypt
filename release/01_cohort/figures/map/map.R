#!/usr/bin/env Rscript
# map.R — where the wheat field cohort was collected.
#
#   Rscript release/01_cohort/scripts/map_cohort.R
# Inputs : map_points.tsv, map_summary.tsv  (prep_map.py)
# Output : map.{png,pdf}
#          map_europe.{png,pdf}
#          map_legend.{png,pdf}
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
# POINT SIZE IS SQRT, so area is proportional to count. Changed 2026-11-26, Leon's call:
# the country marks did not read as a quantity. They were size-coded all along, but log
# put everything from 2 runs upward into the top two-thirds of the range, so 10 runs and
# 65 runs differed by 1.6x in diameter and looked alike. Under sqrt that pair is 2.35x.
#
# The log scale came from the crop map, where counts ran 1 to 1,396 and the small places
# would have vanished. It does not transfer: here the range is 1 to 65, so nothing needs
# protecting at the bottom, and the cost of log is paid at the top where the country
# marks live and where the magnitude is the whole point of drawing them.
#
# What sqrt gives up is the 1-against-2 distinction among localities, which collapses
# from 1.20/2.16 to 1.20/1.54. That is the right thing to spend: most localities are one
# or two runs and the figure says "broad, not deep" by having many small marks, not by
# letting a reader tell one from two.
#
# TWO COLOURS, ONE SHAPE. Green is a named place, blue is a country centroid, for the
# 289 runs whose papers record a country and nothing finer. Diamonds did this job first
# and were swapped out: at these sizes a diamond and a circle are hard to tell apart in
# a cluster, and the difference matters most exactly where the map is densest. Both
# colours are the network figure's, so panels A and B speak the same language.
#
# THE BIGGEST MARKS ARE THE LEAST PRECISE, which is the opposite of what a reader
# expects, so it has to be said rather than designed around. South Africa's 65-run
# circle is the largest symbol on the map and it means "somewhere in South Africa".
# Both series share one size scale so the quantities stay comparable.
#
# A COUNTRY CIRCLE IS NOT A FIELD and the caption has to say so. Those runs were chased hard —
# three papers' supplements read end to end — and the resolution is the literature's,
# not ours. Drawing them as circles would claim a precision that does not exist;
# leaving them off would hide a fifth of the cohort.
#
# THE EUROPEAN WINDOW IS A SECOND IMAGE, NOT AN INSET. Leon's call. It was an inset
# first, and the placement worked, but every property of the box was load-bearing and
# none of it survives a change: the frame proportions had to match the content's own
# coord_fixed or annotation_custom letterboxed it with white bands; the box had to sit
# over ground it itself redrew, so it was widened until it hid no point it did not show
# larger, which pinned it at 72W to 68E; and the frame had to be drawn on the parent
# because annotation_custom eats the outer half of a stroke set on the grob. All of it
# is one figure's worth of tuning that a composition tool does in a drag. Both images
# are written here and assembled on the page.
#
# The two are drawn to ONE SIZE SCALE so the marks stay comparable between them, which
# only holds if they are placed at the same scale on the page. LOCATOR would outline
# the window on the world map; it is OFF, Leon's call, because the frame is the first
# thing the eye lands on and the world map's job is the global spread. Where the second
# image sits is said in the caption instead.
#
# WHAT IS STILL MISSING. 32 runs name no place at all. That is in the summary table.

suppressMessages({
  library(ggplot2); library(maps); library(ragg); library(cowplot)
})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(HERE) || HERE == "") HERE <- "release/01_cohort/figures/map"
pts     <- read.delim(file.path(HERE, "map_points.tsv"),
                      stringsAsFactors = FALSE)
summary <- read.delim(file.path(HERE, "map_summary.tsv"),
                      stringsAsFactors = FALSE)
get <- function(name) summary$value[summary$measure == name]

LEGEND_ON_MAP <- FALSE
LOCATOR       <- FALSE  # outline the European window on the world map

# RENDERED AT FINAL PAGE SIZE, height then width in centimetres, Leon's measurements
# from the assembled page. A ggplot point size is absolute, in millimetres, and does
# NOT scale with the device, so resizing an image on the page changes how big a mark is
# relative to its map. Rendering at the size each one is placed at is the only way the
# key means the same thing on all three, and the only way a 20-run mark is the same
# mark in the world view and the European one.
#
# Changing a number here changes the mark-to-map ratio in that image. It is not a
# resolution setting.
CM        <- 1 / 2.54
MAP_CM    <- c(14.00, 25.67)
EUROPE_CM <- c( 5.65,  8.00)
LEGEND_CM <- c( 7.00,  5.25)

LAND   <- "grey92"      # filled, not outlined: the points are the subject
BORDER <- "white"
POINT   <- "#4a9b34"    # the field-wheat green the network figure uses
COUNTRY <- "#2a78d6"    # the BioProject blue from the same figure
OUTLINE <- "black"

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
  # Breaks reach 65 because the largest marks on the map are country centroids, not
  # localities: a centroid carries 1-65 runs and a named place only 1-30. A key that
  # stopped at 30 left the biggest symbols on the figure unexplained.
  scale_size_continuous(name = "runs", trans = "sqrt", limits = c(1, 65),
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

# EUROPE, NORTH AFRICA AND THE NEAR EAST. 289 of the 460 points and 852 of the runs
# sit in this window, so at world scale most of the cohort is one blot.
#
# The window is cropped to the data, not to the region: south to Tenerife at 28.5N and
# east to Turkey at 35.5E, which are the real corners of what is in it. Rounding out to
# a tidy "Europe and the Mediterranean" box spent a third of the panel on empty Sahara
# and Atlantic and shrank everything worth seeing.
EU <- c(-19, 38, 27, 58)

if (LOCATOR) {
  p <- p + annotate("rect", xmin = EU[1], xmax = EU[2], ymin = EU[3], ymax = EU[4],
                    fill = NA, colour = OUTLINE, linewidth = 0.5)
}

# Same layers and the SAME size scale as the world map, so a mark means the same thing
# in both. It is a second view of one dataset, not a second figure.
europe <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = LAND, colour = BORDER, linewidth = 0.15) +
  geom_point(data = pts, aes(lon, lat, size = runs, fill = placement),
             shape = 21, colour = OUTLINE, stroke = 0.4, alpha = 0.78) +
  scale_fill_manual(values = c("Locality" = POINT, "Country" = COUNTRY),
                    guide = "none") +
  scale_size_continuous(trans = "sqrt", limits = c(1, 65), range = c(1.2, 7),
                        guide = "none") +
  coord_fixed(1.3, xlim = EU[1:2], ylim = EU[3:4], expand = FALSE) +
  theme_void() +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin      = margin(0, 0, 0, 0))

# The canvas has to match the window's own display ratio or the device letterboxes it
# with white bands: 57 degrees of longitude over 31 of latitude at coord_fixed(1.3)
# wants 1.414 to 1, and the placed box is 8.00 / 5.65 = 1.416. That agreement is not a
# coincidence and it is not free; if EUROPE_CM changes, check it again.
stopifnot(abs((EUROPE_CM[2] / EUROPE_CM[1]) /
              (diff(EU[1:2]) / (diff(EU[3:4]) * 1.3)) - 1) < 0.02)

agg_png(file.path(HERE, "map.png"), width = MAP_CM[2] * CM, height = MAP_CM[1] * CM,
        units = "in", res = 600, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "map.pdf"), width = MAP_CM[2] * CM, height = MAP_CM[1] * CM)
print(p); invisible(dev.off())

# The key on its own. Built by drawing the same scales over an empty panel and keeping
# only the guide box, so it cannot drift from what the map actually drew.
key <- ggplot(pts, aes(lon, lat, size = runs, fill = placement)) +
  geom_point(shape = 21, colour = OUTLINE, stroke = 0.4, alpha = 0.78) +
  scale_fill_manual(name = NULL, values = c("Locality" = POINT, "Country" = COUNTRY),
                    labels = c("Locality" = "Locality",
                               "Country" = "Country")) +
  scale_size_continuous(name = "runs", trans = "sqrt", limits = c(1, 65),
                        range = c(1.2, 7), breaks = c(1, 5, 20, 65)) +
  guides(shape = guide_legend(override.aes = list(size = 4, fill = POINT), order = 1),
         size  = guide_legend(override.aes = list(shape = 21, fill = POINT), order = 2)) +
  theme_minimal(base_size = 16) +
  theme(legend.background = element_rect(fill = "white", colour = NA))
# get_plot_component can return a list of candidate grobs; take the first real one.
legend <- cowplot::get_plot_component(key, "guide-box-right", return_all = TRUE)
if (inherits(legend, "list")) legend <- legend[[which(
  vapply(legend, function(g) inherits(g, "gtable") || inherits(g, "grob"), logical(1)))[1]]]

agg_png(file.path(HERE, "map_europe.png"),
        width = EUROPE_CM[2] * CM, height = EUROPE_CM[1] * CM,
        units = "in", res = 600, background = "white")
print(europe); invisible(dev.off())
cairo_pdf(file.path(HERE, "map_europe.pdf"),
          width = EUROPE_CM[2] * CM, height = EUROPE_CM[1] * CM)
print(europe); invisible(dev.off())

agg_png(file.path(HERE, "map_legend.png"),
        width = LEGEND_CM[2] * CM, height = LEGEND_CM[1] * CM,
        units = "in", res = 600, background = "white")
grid::grid.newpage(); grid::grid.draw(legend); invisible(dev.off())
cairo_pdf(file.path(HERE, "map_legend.pdf"),
          width = LEGEND_CM[2] * CM, height = LEGEND_CM[1] * CM)
grid::grid.newpage(); grid::grid.draw(legend); invisible(dev.off())

cat(sprintf("-> map.{png,pdf}  %.2f x %.2f cm; %s points; %s runs at a locality, %s at a country centroid, of %s field-wheat runs\n",
            MAP_CM[1], MAP_CM[2], nrow(pts), get("drawn at a locality"),
            get("drawn at a country centroid"), get("field-wheat runs in cohort")))
cat(sprintf("-> map_europe.{png,pdf}  %.2f x %.2f cm, %s of the %s points inside %s\n",
            EUROPE_CM[1], EUROPE_CM[2],
            sum(pts$lon >= EU[1] & pts$lon <= EU[2] &
                pts$lat >= EU[3] & pts$lat <= EU[4]),
            nrow(pts),
            sprintf("%g..%gE %g..%gN", EU[1], EU[2], EU[3], EU[4])))
cat(sprintf("   %s runs name no place at all and are NOT drawn\n",
            get("not drawn: no usable place")))
