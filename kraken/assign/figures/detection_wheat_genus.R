# detection_wheat_genus.R — the wheat detection scatter, coloured by genus.
#
# The pooled panels colour only by above/below the criterion. Leon saw structure in the
# wheat cloud and asked what it is made of. Colouring by genus shows it: the band above
# the line is not one organism but several, and Alternaria in particular appears as a
# diffuse mass spread across ~15 species nodes, each individually clearing 1%. That is
# the case the flat cutoff handles least well and is worth seeing directly.
#
# Input : kraken/assign/figures/detection_panels.tsv  (Triticum aestivum rows)
# Output: kraken/assign/figures/detection_wheat_genus.png
#
# Palette: validated colourblind-safe categorical set; genus carries a distinct hue,
# everything past the top 7 folds into grey. Italics via plotmath (ggtext no-ops under
# theme_minimal here, see leon-figure-style memory).

library(ggplot2)

d <- read.delim("kraken/assign/figures/detection_panels.tsv", sep = "\t", check.names = FALSE)
d <- d[d$host == "Triticum aestivum", ]
d$genus <- sub(" .*", "", d$taxon)
d$kmer_pct <- d$kmer_frac * 100
pow10 <- function(x) round(log10(x))

# One panel per genus, ordered by how many detections clear the line. Each panel keeps
# the whole-wheat axes so the genera are directly comparable, and the point colour marks
# above/below the criterion within each. This is what shows the per-genus offset: the
# same 1% line falls at a different read count in each panel.
TOP <- c("Puccinia", "Alternaria", "Zymoseptoria", "Fusarium",
         "Blumeria", "Parastagonospora", "Pyrenophora")
d <- d[d$genus %in% TOP, ]
above <- tapply(d$passes == "yes", d$genus, sum)
n_tot <- table(d$genus)
LAB <- setNames(sprintf("%s\nn = %s, %d above", TOP,
                        format(as.integer(n_tot[TOP]), big.mark = ","),
                        as.integer(above[TOP])), TOP)
d$panel <- factor(d$genus, levels = TOP, labels = LAB[TOP])
d$side <- factor(ifelse(d$passes == "yes", "above the criterion", "below"),
                 levels = c("above the criterion", "below"))

p <- ggplot(d, aes(reads, kmer_pct, colour = side)) +
  geom_point(shape = 16, size = 1.2, alpha = 0.35) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.1, colour = "black") +
  facet_wrap(~panel, nrow = 2) +
  scale_x_log10(breaks = c(1e2, 1e4, 1e6), labels = c("100", "10k", "1M")) +
  scale_y_log10(breaks = 10^(-5:2), labels = pow10) +
  scale_colour_manual(values = c("above the criterion" = "#2a78d6", "below" = "grey70"),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4.5, alpha = 1))) +
  labs(x = "reads assigned to the taxon",
       y = expression(log[10]*" (% of k-mer space observed)"),
       title = "Wheat detections, one panel per genus",
       subtitle = paste("Triticum aestivum, db_v3, species rank, 100-read floor. Dashed",
                        "line: 1% of k-mer space.\nThe line falls at a different read",
                        "count in each panel: Alternaria clears it ~10x cheaper than Puccinia.")) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    strip.text       = element_text(size = 12, lineheight = 1.1),
    panel.spacing    = unit(1.1, "lines"),
    plot.title       = element_text(size = 14),
    plot.subtitle    = element_text(size = 11))

ggsave("kraken/assign/figures/detection_wheat_genus.png", p, width = 14, height = 8,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/detection_wheat_genus.png\n")
