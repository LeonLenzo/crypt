#!/usr/bin/env Rscript
# map_localities.R — the cohort at locality resolution, where the metadata allows it.
#
# One point per geocoded locality (OSM Nominatim, see geocode_localities.py): area = samples,
# fill = median collection year. This is the view the country map could not give: the rust
# surveillance corpora resolve to real field sites (PRJEB39201 across Europe, PRJEB31334 across
# Ethiopia) rather than a single national label.
#
# A world panel plus a Europe inset, because 20 of the 53 countries are European and the
# sites there overlap at world scale. Points jitter slightly so coincident sites separate.
#
# Coverage is explicit in the subtitle: 1,450 of 2,643 BioSamples reach a locality; the other
# 794 are country-only (on the country map, not here) and 399 have no location. So this is a
# map of WHERE WE KNOW the sites, not of the whole cohort.
#
#   Rscript 02_literature/03_classify/figures/map_localities.R
# Output: figures/map_localities.{png,pdf}

suppressMessages({
  library(ggplot2); library(maps); library(dplyr); library(scales)
  library(patchwork); library(ragg)
})

here <- "02_literature/03_classify"
pts  <- read.delim(file.path(here, "data/localities_plot.tsv"), stringsAsFactors = FALSE)
prov <- read.delim(file.path(here, "data/cohort_provenance.tsv"), stringsAsFactors = FALSE)

PERIODS <- c("pre-2000", "2000-2009", "2010-2014", "2015-2019", "2020-2025", "unknown")
pts$period <- factor(pts$period, levels = PERIODS)
pts <- pts %>% arrange(desc(n))     # big circles drawn first, small on top

fills <- c("pre-2000" = "#c6dbef", "2000-2009" = "#9ecae1", "2010-2014" = "#6baed6",
           "2015-2019" = "#3182bd", "2020-2025" = "#08519c", "unknown" = "#bdbdbd")

world <- map_data("world"); world <- world[world$region != "Antarctica", ]

base <- theme_minimal(base_size = 15) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

pt_layer <- function(xlim, ylim, maxsize) {
  list(
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "grey94", colour = "grey80", linewidth = 0.25),
    geom_point(data = pts, aes(lon, lat, size = n, fill = period),
               shape = 21, colour = "black", stroke = 0.7, alpha = 0.8,
               position = position_jitter(width = 0.35, height = 0.35, seed = 1)),
    scale_size_area(max_size = maxsize, breaks = c(1, 5, 20, 60),
                    name = "BioSamples"),
    scale_fill_manual(values = fills, name = "Median collection year", drop = TRUE),
    coord_fixed(1.35, xlim = xlim, ylim = ylim, expand = FALSE),
    base,
    theme(axis.title = element_blank(), axis.text = element_blank(),
          panel.grid = element_blank())
  )
}

p_world <- ggplot() + pt_layer(c(-165, 180), c(-56, 78), 14) +
  guides(fill = guide_legend(override.aes = list(size = 5), order = 2),
         size = guide_legend(order = 1)) +
  labs(title = "Cohort sampling localities",
       subtitle = sprintf("%s of %s BioSamples geocoded to %d sites; %s country-only, %s no location",
                          comma(sum(pts$n)), comma(nrow(prov)), nrow(pts),
                          comma(sum(prov$location != "" &
                                    !paste(prov$location, prov$country) %in%
                                      paste(pts$location, pts$country)) -
                                sum(prov$location == "")),
                          comma(sum(prov$location == "")))) +
  theme(legend.position = "right")

p_eu <- ggplot() + pt_layer(c(-11, 32), c(35, 60), 11) +
  guides(fill = "none", size = "none") +
  labs(subtitle = "Europe") +
  theme(plot.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.6))

out <- p_world / p_eu + plot_layout(heights = c(1.5, 1.4))

agg_png(file.path(here, "figures/map_localities.png"), width = 14, height = 13,
        units = "in", res = 300, background = "white")
print(out); invisible(dev.off())
cairo_pdf(file.path(here, "figures/map_localities.pdf"), width = 14, height = 13)
print(out); invisible(dev.off())
cat(sprintf("wrote figures/map_localities.{png,pdf}: %d sites, %d samples\n",
            nrow(pts), sum(pts$n)))
