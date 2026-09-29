# model2d.R — the combined criterion: distinct-minimizer level x accumulation slope.
#
# Neither axis suffices alone. LEVEL (absolute distinct minimizers, = KrakenUniq's unique
# k-mer count) separates the pile-up floor (~1k) from real detections (>10k) regardless of
# read depth. SLOPE separates accumulating from flat. A (host x species) is rejected only
# in the bottom-left: low distinct AND flat. Saturated-real (Phakopsora on soy: low slope,
# huge distinct) is kept by level; real-but-shallow that still accumulates (Bipolaris
# maydis on maize: modest distinct, steep slope) is kept by slope.
#
# Input: kraken/assign/figures/model2d.tsv   Output: model2d.png

library(ggplot2)

d <- read.delim("kraken/assign/figures/model2d.tsv", sep = "\t", check.names = FALSE)
L_CUT <- 1e4; S_CUT <- 0.55
d$verdict <- ifelse(d$med_distinct < L_CUT & d$slope < S_CUT, "reject", "keep")

lab <- c("Phakopsora pachyrhizi|Glycine max|soy Phakopsora (real, saturated)",
         "Phakopsora pachyrhizi|Triticum aestivum|wheat Phakopsora (pile-up)",
         "Bipolaris maydis|Zea mays|B. maydis on maize (real)",
         "Bipolaris maydis|Triticum aestivum|B. maydis on wheat (cross-map)",
         "Puccinia striiformis|Triticum aestivum|stripe rust on wheat (real)",
         "Puccinia striiformis|Zea mays|stripe rust on maize (cross-map)",
         "Sclerotinia sclerotiorum|Brassica napus|Sclerotinia on canola (real)")
lm <- do.call(rbind, lapply(strsplit(lab, "\\|"), function(x)
  data.frame(species = x[1], host = x[2], tag = x[3])))
d <- merge(d, lm, by = c("species", "host"), all.x = TRUE)
key <- d[!is.na(d$tag), ]

p <- ggplot(d, aes(med_distinct, slope)) +
  annotate("rect", xmin = 1, xmax = L_CUT, ymin = -0.3, ymax = S_CUT,
           fill = "#eb6834", alpha = 0.10) +
  annotate("text", x = 30, y = -0.22, hjust = 0, size = 4, colour = "#b3491e",
           label = "reject: low distinct AND flat\n(pile-up / low-read cross-map)") +
  geom_hline(yintercept = S_CUT, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = L_CUT, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = verdict), size = 2, alpha = 0.6) +
  geom_point(data = key, shape = 21, size = 3.4, stroke = 1.1, fill = NA, colour = "black") +
  ggrepel::geom_text_repel(data = key, aes(label = tag), size = 3.3,
                           box.padding = 0.6, min.segment.length = 0, seed = 1,
                           max.overlaps = 20) +
  scale_x_log10(breaks = 10^(2:6), labels = c("100","1k","10k","100k","1M")) +
  scale_colour_manual(values = c(keep = "#2a78d6", reject = "#eb6834"), name = NULL) +
  labs(x = "median distinct minimizers  (level; = KrakenUniq unique k-mers)",
       y = "accumulation slope  (shape)",
       title = "Combined criterion: distinct-minimizer level x accumulation slope",
       subtitle = "Each point a (host x species). Rejected only bottom-left; level rescues saturated reals, slope rescues shallow reals.") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "top", plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 10.5))

ggsave("kraken/assign/figures/model2d.png", p, width = 12, height = 8,
       dpi = 200, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/model2d.png\n")
