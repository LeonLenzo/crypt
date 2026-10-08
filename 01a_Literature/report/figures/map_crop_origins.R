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
# 2,784 runs across 6 localities and wheat is ~1,000 across 364, because the rice data are
# large GWAS panels grown at one farm while the wheat data are rust surveys that sample
# broadly. A map alone therefore reads as "mostly wheat" and a bar chart alone reads as
# "mostly rice". Both panels are present for that reason, and the localities panel exists
# because it is the only one that shows the difference.
#
# Points are drawn smallest-last so a 1,572-run circle cannot hide the single-run localities
# sitting under it. Europe is inset because 2 of every 3 wheat localities fall inside it and
# would otherwise be one blob.
#
# Barley keeps its own colour at three runs. The near-empty key is the finding: there is no
# barley study in the field cereal cohort, only incidental samples inside rust surveys.
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

brks <- c(1, 10, 100, 500)
brks <- brks[brks <= max(pts$n)]
top <- floor(max(pts$n) / 100) * 100
if (top > max(brks)) brks <- c(brks, top)
map_layers <- function(xlim, ylim, ratio, max_size) {
  list(
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "grey94", colour = "grey80", linewidth = 0.3),
    geom_point(data = pts, aes(lon, lat, size = n, fill = crop),
               shape = 21, colour = "black", stroke = 1, alpha = 0.8),
    scale_size_area(max_size = max_size, breaks = brks, labels = comma, name = "Samples"),
    scale_fill_manual(values = fills, name = "Crop", drop = FALSE),
    coord_fixed(ratio, xlim = xlim, ylim = ylim, expand = FALSE)
  )
}

p_map <- ggplot() +
  map_layers(c(-165, 180), c(-48, 72), 1.35, 20) +
  guides(fill = guide_legend(override.aes = list(size = 6), order = 1, ncol = 1),
         size = guide_legend(order = 2)) +
  labs(title = "Field cereal RNA-seq: where the cohort was collected",
       subtitle = sprintf(
         "%s of %s field cereal samples, at %d localities in %d countries. Point area is samples; the remaining %d carry no locality finer than a country.",
         comma(n_pl), comma(n_cer), n_loc, length(unique(pts$country)), n_cer - n_pl)) +
  base +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), legend.position = "right",
        legend.box = "vertical", plot.title = element_text(face = "bold"))

# Europe holds most of the wheat localities and is one blob at world scale. The inset shows
# WHERE they are, at a fixed point size: reusing scale_size_area with a smaller max_size would
# draw the same sample count at two different areas inside one figure, which is a worse lie
# than dropping the encoding.
eur <- pts %>% filter(lon > -12, lon < 42, lat > 34, lat < 62)
p_eur <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey94", colour = "grey80", linewidth = 0.3) +
  geom_point(data = eur, aes(lon, lat, fill = crop), shape = 21, size = 2.6,
             colour = "black", stroke = 0.7, alpha = 0.85) +
  scale_fill_manual(values = fills, drop = FALSE) +
  coord_fixed(1.6, xlim = c(-12, 42), ylim = c(34, 62), expand = FALSE) +
  guides(fill = "none") +
  labs(subtitle = sprintf("Europe: %d of %d localities (points not scaled)",
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
