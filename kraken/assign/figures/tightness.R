# tightness.R — the corrected model: accumulation tightness (R2) vs slope.
#
# Replaces the level x slope model. Level (absolute distinct k-mers) unfairly removed
# low-abundance real detections. Tightness (R2 of log distinct ~ log reads) does not:
# it measures whether coverage tracks depth cleanly, independent of how much there is.
# Real = tight (high R2); pile-up = flat/scattered (low R2); cross-map = diffuse (low R2).
#
# Input: kraken/assign/figures/tightness.tsv   Output: tightness.png

library(ggplot2)
d <- read.delim("kraken/assign/figures/tightness.tsv", sep = "\t", check.names = FALSE)
R_CUT <- 0.5
d$verdict <- ifelse(d$r2 >= R_CUT & d$slope > 0, "keep", "reject")

lab <- c("Puccinia striiformis|Triticum aestivum|stripe rust / wheat (real)",
         "Zymoseptoria tritici|Triticum aestivum|Septoria / wheat (real)",
         "Phakopsora pachyrhizi|Glycine max|soy rust / soybean (real, saturated)",
         "Bipolaris maydis|Zea mays|B. maydis / maize (real, low abundance)",
         "Phytophthora nicotianae|Solanum lycopersicum|P. nicotianae / tomato (real, low abundance)",
         "Phakopsora pachyrhizi|Triticum aestivum|soy rust / wheat (pile-up)",
         "Stemphylium lycopersici|Zea mays|Stemphylium / maize (cross-map)",
         "Bipolaris maydis|Triticum aestivum|B. maydis / wheat (cross-map)")
lm <- do.call(rbind, lapply(strsplit(lab, "\\|"), function(x)
  data.frame(species = x[1], host = x[2], tag = x[3])))
key <- merge(d, lm, by = c("species","host"))

p <- ggplot(d, aes(slope, r2)) +
  annotate("rect", xmin = -0.3, xmax = 1.4, ymin = -0.02, ymax = R_CUT,
           fill = "#eb6834", alpha = 0.09) +
  annotate("text", x = 1.35, y = 0.05, hjust = 1, size = 4, colour = "#b3491e",
           label = "reject: loose fit (pile-up or diffuse cross-map)") +
  geom_hline(yintercept = R_CUT, linetype = "dashed", colour = "grey55") +
  geom_point(aes(colour = verdict), size = 2, alpha = 0.55) +
  geom_point(data = key, shape = 21, size = 3.4, stroke = 1.1, fill = NA, colour = "black") +
  ggrepel::geom_text_repel(data = key, aes(label = tag), size = 3.2,
                           box.padding = 0.7, min.segment.length = 0, seed = 3, max.overlaps = 20) +
  scale_colour_manual(values = c(keep = "#2a78d6", reject = "#eb6834"), name = NULL) +
  labs(x = "accumulation slope  (does coverage grow with depth?)",
       y = "tightness  R²  (does coverage track depth cleanly?)",
       title = "Corrected model: accumulation tightness, not absolute level",
       subtitle = "Each point a (host x species). Tightness is abundance-independent: it keeps low-abundance real\ndetections (Phytophthora/tomato, B. maydis/maize) that a level floor removed, and rejects diffuse cross-maps.") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "top", plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 10.5))

ggsave("kraken/assign/figures/tightness.png", p, width = 12, height = 8,
       dpi = 200, bg = "white", device = ragg::agg_png)
cat("wrote tightness.png\n")
