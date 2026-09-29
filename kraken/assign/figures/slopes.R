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
library(patchwork)

CUT <- 0.55

co <- read.delim("kraken/assign/figures/slopes_cohort.tsv", sep = "\t", check.names = FALSE)
co <- co[order(co$slope), ]
co$species <- factor(co$species, levels = co$species)
co$flag <- ifelse(co$slope < CUT, "artefact-prone", "accumulates")

a <- ggplot(co, aes(slope, species, colour = flag)) +
  geom_vline(xintercept = CUT, linetype = "dashed", linewidth = 1) +
  geom_segment(aes(x = 0, xend = slope, yend = species), linewidth = 0.4, colour = "grey80") +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("artefact-prone" = "#eb6834", "accumulates" = "#2a78d6"),
                      name = NULL) +
  scale_x_continuous(limits = c(0, 1.35), expand = c(0, 0)) +
  labs(x = "accumulation slope   (log distinct minimizers / log reads)", y = NULL,
       title = "A. Pathogen species by accumulation slope (whole cohort)",
       subtitle = sprintf("64 fittable species; %d flagged below %.2f",
                          sum(co$slope < CUT), CUT)) +
  theme_minimal(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.text.y = element_text(size = 7, face = "italic"),
        legend.position = "top", plot.title = element_text(size = 13))

bh <- read.delim("kraken/assign/figures/slopes_byhost.tsv", sep = "\t", check.names = FALSE)
# order species by their spread across hosts (most host-variable at top)
sp_order <- names(sort(tapply(bh$slope, bh$species, function(v) max(v) - min(v))))
bh$species <- factor(bh$species, levels = sp_order)
bh$host <- factor(bh$host, levels = c("Triticum aestivum","Hordeum vulgare","Zea mays",
                  "Oryza sativa","Glycine max","Solanum lycopersicum","Vitis vinifera",
                  "Brassica napus"))

b <- ggplot(bh, aes(host, species, fill = slope)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", slope)), size = 2.7,
            colour = ifelse(bh$slope < CUT, "white", "black")) +
  scale_fill_gradient2(midpoint = CUT, low = "#b3261e", mid = "grey92", high = "#2a78d6",
                       name = "slope", limits = c(0, 1.35)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "B. Slope per host x species",
       subtitle = "blue = accumulates (real), red = plateaus (artefact). A taxon can be real on\nits host and cross-map on another. Midpoint = 0.55 cut.") +
  theme_minimal(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.text.x.top = element_text(angle = 30, hjust = 0, size = 9, face = "italic"),
        axis.text.y = element_text(size = 8, face = "italic"),
        plot.title = element_text(size = 13),
        plot.margin = margin(6, 40, 6, 6))

p <- a + b + plot_layout(widths = c(1, 1.2))
ggsave("kraken/assign/figures/slopes.png", p, width = 17, height = 9,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/slopes.png\n")
