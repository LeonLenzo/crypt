# complexity_byhost.R — reads vs distinct minimizers, one panel per species, points by host.
#
# Companion to the per-host slope heatmap (slopes.png): the heatmap gives the number, this
# shows the cloud behind it. Within a species panel each host is a colour, so a taxon that
# is real on one host and cross-maps on another shows two clouds with different slopes -
# Bipolaris maydis maize points ride the diagonal, wheat points sit on the floor. Per-host
# slope values are printed in each panel. Guide line is distinct/read = 0.03.
#
# Every species that appears in the heatmap gets a panel (37).
# Inputs: complexity_byhost.tsv, complexity_byhost_slopes.tsv   Output: complexity_byhost.png

library(ggplot2)

d  <- read.delim("kraken/assign/figures/complexity_byhost.tsv", sep = "\t", check.names = FALSE)
sl <- read.delim("kraken/assign/figures/complexity_byhost_slopes.tsv", sep = "\t", check.names = FALSE)

HOSTS <- c("Triticum aestivum","Hordeum vulgare","Zea mays","Oryza sativa",
           "Glycine max","Solanum lycopersicum","Vitis vinifera","Brassica napus")
PAL <- setNames(c("#2a78d6","#eb6834","#1baf7a","#eda100",
                  "#e87ba4","#4a3aa7","#008300","#00a0b0"), HOSTS)
d$host  <- factor(d$host, levels = HOSTS)
sl$host <- factor(sl$host, levels = HOSTS)

# panel order: most host-variable species first (widest slope spread across hosts)
spread <- tapply(sl$slope, sl$species, function(v) if (length(v) > 1) max(v) - min(v) else 0)
ord <- names(sort(spread, decreasing = TRUE))
d$species  <- factor(d$species, levels = ord)
sl$species <- factor(sl$species, levels = ord)

# one text line per host-slope, stacked top-left of each panel, coloured by host
sl <- sl[order(sl$species, sl$host), ]
sl$row <- ave(seq_len(nrow(sl)), sl$species, FUN = seq_along)
sl$lab <- sprintf("%.2f", sl$slope)

p <- ggplot(d, aes(reads, distinct, colour = host)) +
  geom_abline(slope = 1, intercept = log10(0.03), linetype = "dashed",
              linewidth = 0.7, colour = "grey40") +
  geom_point(size = 0.7, alpha = 0.35) +
  geom_text(data = sl, aes(x = 100, y = 10^(7 - 0.55 * (row - 1)), label = lab, colour = host),
            hjust = 0, size = 3, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~species, ncol = 6) +
  scale_x_log10(breaks = 10^c(2,4,6), labels = c("100","10k","1M")) +
  scale_y_log10(breaks = 10^c(1,4,7), labels = c("10","10k","10M")) +
  scale_colour_manual(values = PAL, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3.5, alpha = 1))) +
  labs(x = "reads assigned to the taxon", y = "distinct minimizers observed",
       title = "Reads vs distinct minimizers per species, coloured by host",
       subtitle = paste("Numbers = per-host accumulation slope. Dashed = distinct/read =",
                        "0.03. A species real on its host and cross-mapping on another\nshows",
                        "clouds with different slopes (e.g. Bipolaris maydis: maize diagonal,",
                        "wheat floor).")) +
  theme_minimal(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(size = 9, face = "italic"),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "top", legend.text = element_text(size = 10),
        plot.title = element_text(size = 15), plot.subtitle = element_text(size = 11))

ggsave("kraken/assign/figures/complexity_byhost.png", p, width = 16, height = 15,
       dpi = 200, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/complexity_byhost.png\n")
