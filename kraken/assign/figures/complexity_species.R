# complexity_species.R — reads vs distinct minimizers, one panel per species.
#
# The pooled genus scatter mixed species with different genome sizes, so the clouds
# overlapped. Per species the relationship is clean: a genuinely present organism rides
# the diagonal (distinct k-mers grow with reads), a pile-up sits on a flat floor
# (distinct pinned while reads climb). The dashed guide is distinct/read = 0.03, the
# bimodal valley from complexity.png. Denominator-free, so genome size does not distort
# the comparison across panels.
#
# Top 20 wheat pathogen species by detection count.
# Input: kraken/assign/figures/complexity_species.tsv   Output: complexity_species.png

library(ggplot2)

d <- read.delim("kraken/assign/figures/complexity_species.tsv", sep = "\t", check.names = FALSE)
THRESH <- 0.03

# order panels by detection count, and count what clears the guide per species
ord <- names(sort(table(d$taxon), decreasing = TRUE))
above <- tapply(d$distinct / d$reads >= THRESH, d$taxon, sum)
n_tot <- table(d$taxon)
lab <- setNames(sprintf("%s\n%d of %d above", ord,
                        as.integer(above[ord]), as.integer(n_tot[ord])), ord)
d$panel <- factor(d$taxon, levels = ord, labels = lab[ord])
d$side <- ifelse(d$distinct / d$reads >= THRESH, "above", "below")

p <- ggplot(d, aes(reads, distinct, colour = side)) +
  geom_abline(slope = 1, intercept = log10(THRESH),
              linetype = "dashed", linewidth = 0.9, colour = "black") +
  geom_point(size = 0.9, alpha = 0.4) +
  facet_wrap(~panel, ncol = 5) +
  scale_x_log10(breaks = 10^c(2,4,6), labels = c("100","10k","1M")) +
  scale_y_log10(breaks = 10^c(1,3,5,7), labels = c("10","1k","100k","10M")) +
  scale_colour_manual(values = c(above = "#2a78d6", below = "grey70"),
                      name = NULL, labels = c("above distinct/read = 0.03", "below")) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(x = "reads assigned to the taxon", y = "distinct minimizers observed",
       title = "Reads vs distinct minimizers, top 20 wheat pathogen species",
       subtitle = paste("Triticum aestivum, db_v3, 100-read floor. Diagonal =",
                        "real (k-mers accumulate); flat floor = pile-up.",
                        "Dashed: distinct/read = 0.03.")) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(size = 9.5, lineheight = 1.05),
        panel.spacing = unit(0.7, "lines"),
        legend.position = "top",
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 11))

ggsave("kraken/assign/figures/complexity_species.png", p, width = 15, height = 8.5,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/complexity_species.png\n")
