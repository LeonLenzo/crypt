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

TOP <- c("Puccinia", "Alternaria", "Zymoseptoria", "Fusarium",
         "Blumeria", "Parastagonospora", "Pyrenophora")
PAL <- c(Puccinia = "#2a78d6", Alternaria = "#eb6834", Zymoseptoria = "#1baf7a",
         Fusarium = "#eda100", Blumeria = "#e87ba4", Parastagonospora = "#4a3aa7",
         Pyrenophora = "#111111", Other = "grey78")
d$grp <- ifelse(d$genus %in% TOP, d$genus, "Other")
d$grp <- factor(d$grp, levels = c(TOP, "Other"))
# draw Other first (underneath), coloured genera on top
d <- d[order(d$grp == "Other", decreasing = TRUE), ]

lab <- setNames(lapply(levels(d$grp),
                       function(g) if (g == "Other") "Other"
                                   else bquote(italic(.(g)))), levels(d$grp))

p <- ggplot(d, aes(reads, kmer_pct, fill = grp)) +
  geom_point(shape = 21, size = 1.9, stroke = 0.15, colour = "grey15", alpha = 0.55) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2, colour = "black") +
  scale_x_log10(breaks = c(1e2, 1e4, 1e6), labels = c("100", "10k", "1M")) +
  scale_y_log10(breaks = 10^(-5:2), labels = pow10) +
  scale_fill_manual(values = PAL, breaks = levels(d$grp),
                    labels = lab, name = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             override.aes = list(size = 4.5, alpha = 1))) +
  labs(x = "reads assigned to the taxon",
       y = expression(log[10]*" (% of k-mer space observed)"),
       title = "Wheat detections by genus",
       subtitle = paste("Triticum aestivum, db_v3, species rank, 100-read floor.",
                        "Dashed line: 1% of k-mer space.\nAlternaria is one genus",
                        "spread across ~15 species nodes, each clearing the line alone.")) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    plot.title       = element_text(size = 14),
    plot.subtitle    = element_text(size = 11))

ggsave("kraken/assign/figures/detection_wheat_genus.png", p, width = 11, height = 7.5,
       dpi = 300, bg = "white", device = ragg::agg_png)
cat("wrote kraken/assign/figures/detection_wheat_genus.png\n")
