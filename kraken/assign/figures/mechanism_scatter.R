# mechanism_scatter.R — why the k-mer criterion separates real detections from artefacts.
#
# One point per (run, taxon) detection passing a 100-read floor: reads on x, the
# fraction of that taxon's k-mer space observed on y, both log. The two axes disagree,
# and that disagreement is the whole argument. Phakopsora sits far right and far down,
# hundreds of thousands of reads landing on a sliver of its genome, while Puccinia
# triticina sits left and high on far fewer reads spread broadly. A vertical read
# threshold cannot separate them; a horizontal one can.
#
# Input : kraken/assign/figures/mechanism_scatter.tsv  (mech.py; the unhighlighted
#         cloud is thinned to 25% deliberately, it is background)
# Output: kraken/assign/figures/mechanism_scatter.{png,pdf}
#
# Self-contained mode. Italics via plotmath, not ggtext: see artefact_collapse.R.

library(ggplot2)

d <- read.delim("kraken/assign/figures/mechanism_scatter.tsv", sep = "\t", check.names = FALSE)

# percent, so the decades land on 10^-3..10^2 and label as powers rather than as
# strings of leading zeros
d$kmer_pct <- d$kmer_frac * 100
# bare exponents, not 10^x: the axis title carries the log10, so the ticks only
# need to say which decade
pow10 <- function(x) round(log10(x))

# NOT "artefact": of the 124 Phakopsora/Melampsora detections above the line, 123 are
# the declaring study's own pathogen. The taxon is not spurious, the low-coverage
# detections of it are, and the criterion tells them apart without seeing the metadata.
LAB <- c(target    = "P. striiformis (declared target)",
         secondary = "P. triticina (genuine secondary)",
         artefact  = "Phakopsora / Melampsora",
         other     = "all other taxa")
# Colourblind-safe categorical palette, validated (lightness band, chroma floor,
# CVD separation, normal-vision floor). NOT the Lenzo et al. 2026 field-data palette:
# those hues are specific to that manuscript's organisms and do not carry over.
PAL <- c(target    = "#2a78d6",
         secondary = "#eb6834",
         artefact  = "#1baf7a",
         other     = "grey75")
SHP <- c(target = 21, secondary = 22, artefact = 24, other = 21)

d$class <- factor(d$class, levels = c("target", "secondary", "artefact", "other"))
d <- d[order(d$class, decreasing = TRUE), ]      # draw the grey cloud underneath

p <- ggplot(d, aes(x = reads, y = kmer_pct)) +
  geom_point(aes(fill = class, shape = class, size = class, alpha = class),
             colour = "black", stroke = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2, colour = "black") +
  annotate("text", x = 2.5e7, y = 2.8, hjust = 1, size = 4.4,
           label = "1% of k-mer space") +
  scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7),
                labels = c("100", "1k", "10k", "100k", "1M", "10M")) +
  # data runs to log10 = -4.74, so the ticks have to reach -5 or the axis is
  # labelled over only part of its range
  scale_y_log10(breaks = 10^(-5:2), labels = pow10) +
  scale_fill_manual(values = PAL, labels = LAB, name = NULL) +
  scale_shape_manual(values = SHP, labels = LAB, name = NULL) +
  scale_size_manual(values = c(target = 2.6, secondary = 2.6, artefact = 2.6, other = 1.6),
                    labels = LAB, name = NULL, guide = "none") +
  scale_alpha_manual(values = c(target = 0.7, secondary = 0.7, artefact = 0.7, other = 0.25),
                     labels = LAB, name = NULL, guide = "none") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             override.aes = list(size = 4, alpha = 0.9)),
         shape = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "reads assigned to the taxon", y = "fraction of the taxon's k-mer space observed",
       title = "Genome coverage separates real detections from read pile-ups",
       subtitle = paste0("3,220 runs across 73 host species (46% wheat), db_v3, species rank;\nother taxa thinned to 25% for legibility. ",
                         "Of 124 Phakopsora/Melampsora detections above the line, 123 are the ",
                         "declaring study's own pathogen.")) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    plot.title       = element_text(size = 13),
    plot.subtitle    = element_text(size = 11))

ggsave("kraken/assign/figures/mechanism_scatter.png", p, width = 11, height = 7.5,
       dpi = 300, bg = "white", device = ragg::agg_png)
ggsave("kraken/assign/figures/mechanism_scatter.pdf", p, width = 11, height = 7.5,
       bg = "white", device = cairo_pdf)
cat("wrote kraken/assign/figures/mechanism_scatter.{png,pdf}\n")
