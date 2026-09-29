# artefact_collapse.R — what the k-mer criterion removes, and from which taxa.
#
# Each taxon gets two points: the runs it was detected in under a read-count floor
# alone (>= 100 clade reads), and the runs surviving the k-mer criterion (>= 1% of
# that taxon's minimizers in db_v3 observed). Both computed over the same 3,220
# reports, so the pair is comparable.
#
# The point of the figure is that the collapse is not arbitrary. Taxa lose runs in
# proportion to how little of their genome the reads ever covered, so the organisms
# that fall away are the ones that were never really there.
#
# Input : kraken/assign/figures/artefact_collapse.tsv  (collapse.py)
# Output: kraken/assign/figures/artefact_collapse.png
#
# Self-contained mode: axis titles, labels and legend kept.

library(ggplot2)

TOP_N <- 16

# base read.delim, not readr: this machine's vroom.so fails to load against R 4.6
d <- read.delim("kraken/assign/figures/artefact_collapse.tsv", sep = "\t", check.names = FALSE)
d <- head(d[order(-d$runs_readfloor), ], TOP_N)

# Italic binomials via plotmath, NOT ggtext. element_markdown() renders literally
# whenever a complete theme is in play, and theme_minimal is the house theme, so it is
# not avoidable here. Upgrading ggtext 0.1.2 -> 0.2.0 does not help: the element really
# is class element_markdown, ggtext's draw method just is not dispatched. No warning.
d$label <- factor(d$taxon, levels = rev(d$taxon[order(d$retained_pct)]))
italicise <- function(x) parse(text = sprintf('italic("%s")', x))

long <- rbind(
  data.frame(label = d$label, runs = d$runs_readfloor, stage = "Read-count floor only (≥ 100 reads)"),
  data.frame(label = d$label, runs = d$runs_kmer,      stage = "k-mer criterion (≥ 1% of k-mer space)")
)
long$stage <- factor(long$stage, levels = c("Read-count floor only (≥ 100 reads)",
                                            "k-mer criterion (≥ 1% of k-mer space)"))

p <- ggplot() +
  geom_segment(data = d,
               aes(y = label, yend = label, x = runs_kmer, xend = runs_readfloor),
               colour = "black", linewidth = 1.2) +
  geom_point(data = long, aes(y = label, x = runs, fill = stage),
             shape = 21, size = 4, stroke = 1, colour = "black") +
  # validated categorical blue, not the Lenzo et al. 2026 manuscript palette
  scale_fill_manual(values = c("Read-count floor only (≥ 100 reads)"    = "grey80",
                               "k-mer criterion (≥ 1% of k-mer space)" = "#2a78d6")) +
  # Median k-mer fraction as a right-hand axis, not a geom_text at an invented x
  # position: padding the x scale to make room would imply runs that do not exist.
  # Putting the mechanism beside the collapse is the whole point of the figure.
  scale_y_discrete(labels = italicise,
                   sec.axis = dup_axis(labels = sprintf("%.3f%%", d$median_kmer_frac[
                     match(levels(d$label), d$taxon)]),
                     name = "median fraction of k-mer space observed")) +
  labs(x = "SRA runs with the taxon detected", y = NULL, fill = NULL,
       title = "Detections before and after the k-mer criterion",
       subtitle = "3,220 runs across 73 host species (46% wheat), db_v3, species rank") +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.position  = "top",
    plot.title       = element_text(size = 13),
    plot.subtitle    = element_text(size = 11),
    axis.text.y.right = element_text(size = 13),
    axis.title.y.right = element_text(size = 12, angle = 90))

# ragg, not the default png device: the default substitutes ">=" for the U+2265 glyph
ggsave("kraken/assign/figures/artefact_collapse.png", p, width = 11, height = 7, dpi = 300,
       bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/artefact_collapse.png\n")
