# complexity.R — the transform that separates real detections from pile-ups.
#
# Left : reads vs distinct minimizers, log-log. A real detection accumulates k-mers as
#        reads grow, so it rides a diagonal; a pile-up saturates, so it lies on a flat
#        floor. The diagonal guide is distinct/read = 0.03 (the histogram valley on the
#        right). This is library complexity, and it needs no genome-size denominator.
# Right: the distribution of log10(distinct/read). It is bimodal, and the valley near
#        -1.5 is an objective threshold rather than a chosen one.
#
# Points coloured by whether the taxon is a known wheat-context artefact, to show the
# separation is real, not assumed.
#
# Input : kraken/assign/figures/complexity.tsv   Output: complexity.png

library(ggplot2)
library(patchwork)

d <- read.delim("kraken/assign/figures/complexity.tsv", sep = "\t", check.names = FALSE)
d$logratio <- log10(d$ratio)
THRESH <- 0.03                      # ~ histogram valley, distinct per read
PAL <- c("other wheat pathogens" = "#2a78d6",
         "artefact taxa (wheat context)" = "#eb6834")

a <- ggplot(d, aes(reads, distinct, colour = kind)) +
  geom_point(size = 1.1, alpha = 0.3) +
  geom_abline(slope = 1, intercept = log10(THRESH), linetype = "dashed",
              linewidth = 1.1, colour = "black") +
  annotate("text", x = 3e6, y = 3e6 * THRESH, label = "distinct/read = 0.03",
           hjust = 1, vjust = -0.6, size = 4) +
  scale_x_log10(breaks = 10^(2:7), labels = c("100","1k","10k","100k","1M","10M")) +
  scale_y_log10(breaks = 10^(0:7),
                labels = c("1","10","100","1k","10k","100k","1M","10M")) +
  scale_colour_manual(values = PAL, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(x = "reads assigned to the taxon", y = "distinct minimizers observed",
       title = "Real detections accumulate k-mers; pile-ups do not") +
  theme_minimal(base_size = 15) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "top", plot.title = element_text(size = 13))

b <- ggplot(d, aes(logratio, fill = kind)) +
  geom_histogram(bins = 40, colour = "grey30", linewidth = 0.15, alpha = 0.85,
                 position = "stack") +
  geom_vline(xintercept = log10(THRESH), linetype = "dashed", linewidth = 1.1) +
  scale_fill_manual(values = PAL, name = NULL, guide = "none") +
  labs(x = expression(log[10]*" (distinct minimizers per read)"), y = "detections",
       title = "Bimodal: the valley is the threshold") +
  theme_minimal(base_size = 15) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(size = 13))

p <- a + b + plot_layout(widths = c(1.4, 1))
ggsave("kraken/assign/figures/complexity.png", p, width = 15, height = 6.5,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/complexity.png\n")
