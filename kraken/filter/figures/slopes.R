# slopes.R — the accumulation-slope filter, in two views.
#
# A: every fittable pathogen species ranked by cohort slope, with the 0.55 cut. Below
#    the cut = plateaus across the cohort = artefact-prone. The flagged set is off-target
#    / non-crop / cross-mapping, found without being told what a crop cohort holds.
# B: slope per (host x species) as a heatmap, for species seen on >=2 major hosts. This
#    is why the filter MUST be per-host: a taxon can be real on its own host (high, blue)
#    and cross-map on another (low, red). Bipolaris maydis is the clean case.
#
# Inputs: slopes_cohort.tsv, slopes_byhost.tsv   Output: slopes.png

library(ggplot2)

CUT <- 0.55

bh <- read.delim("kraken/filter/data/model2d.tsv", sep = "\t", check.names = FALSE)
# cells the combined rule keeps despite a low slope, because coverage (level) is high:
# saturated real detections. These would read as "artefact" on slope colour alone.
L_CUT <- 1e4
bh$rescued <- bh$slope < CUT & bh$med_distinct >= L_CUT
# order species by their spread across hosts (most host-variable at top)
sp_order <- names(sort(tapply(bh$slope, bh$species, function(v) max(v) - min(v))))
bh$species <- factor(bh$species, levels = sp_order)
bh$host <- factor(bh$host, levels = c("Triticum aestivum","Hordeum vulgare","Zea mays",
                  "Oryza sativa","Glycine max","Solanum lycopersicum","Vitis vinifera",
                  "Brassica napus"))

b <- ggplot(bh, aes(host, species, fill = slope)) +
  geom_tile(colour = "grey96", linewidth = 0.6) +
  geom_point(data = function(x) x[x$rescued, ], shape = 21, size = 9, stroke = 1.4,
             fill = NA, colour = "grey15") +
  geom_text(aes(label = sprintf("%.2f", slope),
                colour = abs(slope - CUT) > 0.28), size = 2.8, fontface = "bold") +
  scale_fill_gradient2(midpoint = CUT, low = "#c0392b", mid = "#f4f4f2", high = "#2a78d6",
                       name = "slope", limits = c(-0.15, 1.05),
                       breaks = c(0, CUT, 1), labels = c("0.0", "0.55", "1.0")) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Accumulation slope per host x pathogen species",
       subtitle = "Fill = accumulation slope (blue real, red artefact). Ringed cells are kept anyway:\nlow slope but high genome coverage = a real, fully-covered infection. Cut = 0.55.") +
  theme_minimal(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, face = "italic"),
        axis.text.y = element_text(size = 8, face = "italic"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 13),
        plot.margin = margin(6, 20, 6, 6))

ggsave("kraken/filter/figures/slopes.png", b, width = 10.5, height = 9,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("wrote kraken/filter/figures/slopes.png\n")
