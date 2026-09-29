# detection_panels.R — what a known-real detection looks like, and how the rest compare.
#
# Two panels, split on whether the declaring study named that pathogen, not on species.
# The left panel is a positive control: organisms the authors said were present. The
# right panel is everything else, and it is where the cryptic co-infections live.
#
# The point of the split is that the criterion was never shown the metadata. It keeps
# 94% of what studies declare and rejects 89% of what they do not, so the 2,550
# undeclared detections that survive are candidates rather than leftovers.
#
# Input : kraken/assign/figures/detection_panels.tsv  (prep_detection_panels.py)
# Output: kraken/assign/figures/detection_panels.{png,pdf}
#
# Self-contained mode. Palette is the validated colourblind-safe categorical set.

library(ggplot2)

d <- read.delim("kraken/assign/figures/detection_panels.tsv", sep = "\t", check.names = FALSE)

# Work in percent so the decade labels land on 10^-3 to 10^2, and label them as powers
# rather than as 0.001%/0.1%/1%: on a log axis the exponent is the quantity being
# varied, and reading four leading zeros is harder than reading -3.
d$kmer_pct <- d$kmer_frac * 100
pow10 <- function(x) parse(text = sprintf("10^%d", round(log10(x))))

ORD <- c("declared", "secondary pathogen", "secondary non-pathogen")
TITLE <- c(declared                 = "Declared by the study",
           `secondary pathogen`     = "Secondary, a known pathogen",
           `secondary non-pathogen` = "Secondary, not a known pathogen")
n   <- table(d$class)
pct <- tapply(d$passes == "yes", d$class, mean) * 100
LAB <- sapply(ORD, function(k)
  sprintf("%s\nn = %s, %.1f%% above the criterion",
          TITLE[[k]], format(n[[k]], big.mark = ","), pct[[k]]))
d$class <- factor(d$class, levels = ORD, labels = LAB)
d$passes <- factor(d$passes, levels = c("yes", "no"),
                   labels = c("above the criterion", "below"))

p <- ggplot(d, aes(x = reads, y = kmer_pct)) +
  # No outline here, against the usual house mark. At 52k points per panel a border
  # fuses into a solid slab and the panel stops reading as a density; unoutlined
  # semi-transparent points let the shape of the cloud through instead.
  geom_point(aes(colour = passes), shape = 16, size = 1.5, alpha = 0.45) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2, colour = "black") +
  facet_wrap(~class, nrow = 1) +
  scale_x_log10(breaks = c(1e2, 1e4, 1e6), labels = c("100", "10k", "1M")) +
  scale_y_log10(breaks = 10^(-3:2), labels = pow10) +
  scale_colour_manual(values = c("above the criterion" = "#2a78d6", "below" = "grey70"),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4.5, alpha = 1))) +
  labs(x = "reads assigned to the taxon",
       y = "% of the taxon's k-mer space observed",
       title = "The criterion agrees with the studies without being shown them",
       subtitle = paste0("3,220 runs, db_v3, species rank, 100-read floor. Dashed line: 1% of ",
                         "k-mer space. Host taxa excluded.\nSecondary pathogens and ",
                         "non-pathogens behave alike (11.4% vs 11.3%): the criterion tests ",
                         "presence, not whether an organism is on a pathogen list.")) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(size = 13, lineheight = 1.15),
    legend.position  = "top",
    plot.title       = element_text(size = 14),
    plot.subtitle    = element_text(size = 11),
    panel.spacing    = unit(1.4, "lines"),
    plot.margin      = margin(8, 12, 8, 14))

ggsave("kraken/assign/figures/detection_panels.png", p, width = 14.5, height = 6.2,
       dpi = 300, bg = "white", device = ragg::agg_png)
ggsave("kraken/assign/figures/detection_panels.pdf", p, width = 14.5, height = 6.2,
       bg = "white", device = cairo_pdf)
cat("wrote kraken/assign/figures/detection_panels.{png,pdf}\n")
