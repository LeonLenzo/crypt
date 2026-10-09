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
# Points are NOT scaled by sample count (leon, 2026-10-08: "the bars already do that"). One
# dot is one locality, so the map answers WHERE and the bars answer HOW MANY, and neither
# tries to do both. Scaling also made three 500-to-1,400-sample circles swamp the 378
# single-digit wheat localities that are the point of the panel.
#
# Europe is inset because 2 of every 3 wheat localities fall inside it and would otherwise be
# one blob.
#
# Crops under 100 samples are dropped by prep_crop_map.py, not here, so the bar panels and the
# map always agree about what exists. What was dropped is in crop_totals.tsv with an
# `included` flag; as of 2026-10-08 that is barley (3) and the small cereals (53).
#
#   Rscript 01a_Literature/report/figures/map_crop_origins.R
# Inputs : data/gold/crop_map.tsv, data/gold/crop_totals.tsv  (prep_crop_map.py)
# Output : report/figures/map_crop_origins.{png,pdf}

suppressMessages({
  library(ggplot2); library(maps); library(dplyr); library(scales)
  library(patchwork); library(ragg)
})

here <- "01a_Literature"
pts  <- read.delim(file.path(here, "data/gold/crop_map.tsv"), stringsAsFactors = FALSE)
tot  <- read.delim(file.path(here, "data/gold/crop_totals.tsv"), stringsAsFactors = FALSE)

tot <- tot[tot$included == "yes", ]
ORDER <- c("Wheat", "Maize", "Rice", "Sorghum", "Barley", "Other cereal",
           "Cereal, host unresolved")
ORDER <- ORDER[ORDER %in% unique(c(pts$crop, tot$crop))]
pts$crop <- factor(pts$crop, levels = ORDER)
tot$crop <- factor(tot$crop, levels = ORDER)

# Okabe-Ito, which is colourblind-safe by construction, assigned so that the four crops
# carrying almost all the data are maximally separated. Grey sits outside the ramp for the
# runs whose host the archive never recorded.
fills <- c("Wheat"                   = "#0072B2",
           "Maize"                   = "#E69F00",
           "Rice"                    = "#009E73",
           "Sorghum"                 = "#D55E00",
           "Barley"                  = "#CC79A7",
           "Other cereal"            = "#56B4E9",
           "Cereal, host unresolved" = "grey60")

base <- theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

world <- map_data("world")
world <- world[world$region != "Antarctica", ]

# smallest drawn last, so single-run localities stay visible over the big panels
pts <- pts %>% arrange(desc(n))

# Round the top break DOWN: ggplot drops any break outside the data range, so signif() on
# the maximum (1,396 -> 1,400) silently removed the largest key and left the biggest circle
# on the map with nothing in the legend to measure it against.
n_cer <- sum(tot$runs); n_loc <- length(unique(pts$location)); n_pl <- sum(pts$n)

map_layers <- function(xlim, ylim, ratio, psize) {
  list(
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "grey94", colour = "grey80", linewidth = 0.3),
    geom_point(data = pts, aes(lon, lat, fill = crop), shape = 21, size = psize,
               colour = "black", stroke = 0.8, alpha = 0.85),
    scale_fill_manual(values = fills, name = "Crop", drop = FALSE),
    coord_fixed(ratio, xlim = xlim, ylim = ylim, expand = FALSE)
  )
}

p_map <- ggplot() +
  map_layers(c(-165, 180), c(-48, 72), 1.35, 3.1) +
  guides(fill = guide_legend(override.aes = list(size = 6), ncol = 1)) +
  labs(title = "Field cereal RNA-seq: where the cohort was collected",
       subtitle = sprintf(
         "%s of %s field cereal samples, at %d localities in %d countries. One dot is one locality; sample counts are in the bars below.",
         comma(n_pl), comma(n_cer), n_loc, length(unique(pts$country)))) +
  base +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), legend.position = "right",
        legend.box = "vertical", plot.title = element_text(face = "bold"))

# Europe holds most of the wheat localities and is one blob at world scale, so it gets an
# inset at the same fixed point size as the main map. Nothing in this figure encodes sample
# count as area any more; the bars do that.
eur <- pts %>% filter(lon > -12, lon < 42, lat > 34, lat < 62)
p_eur <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey94", colour = "grey80", linewidth = 0.3) +
  geom_point(data = eur, aes(lon, lat, fill = crop), shape = 21, size = 2.6,
             colour = "black", stroke = 0.7, alpha = 0.85) +
  scale_fill_manual(values = fills, drop = FALSE) +
  coord_fixed(1.6, xlim = c(-12, 42), ylim = c(34, 62), expand = FALSE) +
  guides(fill = "none") +
  labs(subtitle = sprintf("Europe: %d of the %d localities",
                          length(unique(eur$location)), n_loc)) +
  base +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 1))

p_runs <- ggplot(tot, aes(runs, crop, fill = crop)) +
  geom_col(colour = "black", linewidth = 0.9, width = 0.74) +
  geom_text(aes(label = comma(runs)), hjust = -0.2, size = 4.4, colour = "grey20") +
  scale_fill_manual(values = fills, guide = "none") +
  scale_x_continuous(expand = expansion(c(0, 0.22)), labels = comma) +
  scale_y_discrete(limits = rev(ORDER)) +
  labs(x = "Samples", y = NULL, subtitle = "Samples per crop") +
  base + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

p_locs <- ggplot(tot, aes(localities, crop, fill = crop)) +
  geom_col(colour = "black", linewidth = 0.9, width = 0.74) +
  geom_text(aes(label = comma(localities)), hjust = -0.2, size = 4.4, colour = "grey20") +
  scale_fill_manual(values = fills, guide = "none") +
  scale_x_continuous(expand = expansion(c(0, 0.22)), labels = comma) +
  scale_y_discrete(limits = rev(ORDER)) +
  labs(x = "Distinct localities", y = NULL,
       subtitle = "Localities per crop — the same data, reshaped") +
  base + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
               axis.text.y = element_blank())

out <- p_map / (p_eur | p_runs | p_locs) +
  plot_layout(heights = c(2.1, 1.25), widths = c(1, 1, 1))

agg_png(file.path(here, "report/figures/map_crop_origins.png"), width = 17, height = 12,
        units = "in", res = 300, background = "white")
print(out); invisible(dev.off())
cairo_pdf(file.path(here, "report/figures/map_crop_origins.pdf"), width = 17, height = 12)
print(out); invisible(dev.off())
cat("wrote report/figures/map_crop_origins.{png,pdf}\n")
cat(sprintf("  %d locality x crop points, %s samples, %d localities\n",
            nrow(pts), comma(n_pl), n_loc))
