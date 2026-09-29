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

ALL <- read.delim("kraken/assign/figures/detection_panels.tsv", sep = "\t", check.names = FALSE)

# One image per host with enough samples to be worth reading, plus one pooled. The
# cohort spans 73 host species, and a pathogen that is ordinary on one host reads as
# noise when pooled with another's.
MIN_DETECTIONS <- 400
tab   <- sort(table(ALL$host), decreasing = TRUE)
hosts <- names(tab)[tab >= MIN_DETECTIONS & names(tab) != "unresolved"]
jobs  <- c(list(list(key = "all", host = NULL)),
           lapply(hosts, function(h) list(key = gsub("[^A-Za-z0-9]+", "_", tolower(h)), host = h)))

for (job in jobs) {
d <- if (is.null(job$host)) ALL else ALL[ALL$host == job$host, ]

# Work in percent so the decade labels land on 10^-3 to 10^2, and label them as powers
# rather than as 0.001%/0.1%/1%: on a log axis the exponent is the quantity being
# varied, and reading four leading zeros is harder than reading -3.
d$kmer_pct <- d$kmer_frac * 100
# bare exponents, not 10^x: the axis title carries the log10, so the ticks only
# need to say which decade
pow10 <- function(x) round(log10(x))

ORD <- c("declared", "secondary pathogen", "secondary non-pathogen")
TITLE <- c(declared                 = "Declared by the study",
           `secondary pathogen`     = "Secondary, a known pathogen",
           `secondary non-pathogen` = "Secondary, not a known pathogen")
# fixed levels, because a single host can have zero detections in a class and an
# absent level would otherwise drop out of the table and break the labels
n   <- table(factor(d$class, levels = ORD))
pct <- sapply(ORD, function(k) {
  v <- d$passes[d$class == k]
  if (length(v) == 0) NA_real_ else mean(v == "yes") * 100
})
LAB <- sapply(ORD, function(k)
  if (is.na(pct[[k]]))
    sprintf("%s\nno detections", TITLE[[k]])
  else
    sprintf("%s\nn = %s, %.1f%% above the criterion",
            TITLE[[k]], format(n[[k]], big.mark = ","), pct[[k]]))
d$class <- factor(d$class, levels = ORD, labels = LAB)
d <- d[!is.na(d$class), ]
d$passes <- factor(d$passes, levels = c("yes", "no"),
                   labels = c("above the criterion", "below"))

p <- ggplot(d, aes(x = reads, y = kmer_pct)) +
  # No outline here, against the usual house mark. At 52k points per panel a border
  # fuses into a solid slab and the panel stops reading as a density; unoutlined
  # semi-transparent points let the shape of the cloud through instead.
  geom_point(aes(colour = passes), shape = 16, size = 1.4, alpha = 0.18) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2, colour = "black") +
  facet_wrap(~class, nrow = 1) +
  scale_x_log10(breaks = c(1e2, 1e4, 1e6), labels = c("100", "10k", "1M")) +
  # data runs to log10 = -4.74, so the ticks have to reach -5 or the axis is
  # labelled over only part of its range
  scale_y_log10(breaks = 10^(-5:2), labels = pow10) +
  scale_colour_manual(values = c("above the criterion" = "#2a78d6", "below" = "grey70"),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4.5, alpha = 1))) +
  labs(x = "reads assigned to the taxon",
       y = expression(log[10]*" (% of k-mer space observed)"),
       title = if (is.null(job$host)) "The criterion agrees with the studies without being shown them"
               else sprintf("Detections in %s", job$host),
       subtitle = if (is.null(job$host))
                    paste0("3,220 runs across 73 host species (46% wheat), db_v3, species rank, ",
                           "100-read floor. Dashed line: 1% of k-mer space. Host taxa excluded.")
                  else
                    sprintf(paste0("db_v3, species rank, 100-read floor, host taxa excluded. ",
                                   "Dashed line: 1%% of k-mer space.\n%s detections from this host."),
                            format(nrow(d), big.mark = ","))) +
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

stem <- if (is.null(job$host)) "kraken/assign/figures/detection_panels"
        else sprintf("kraken/assign/figures/detection_panels_%s", job$key)
ggsave(paste0(stem, ".png"), p, width = 14.5, height = 6.2,
       dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem, ".pdf"), p, width = 14.5, height = 6.2,
       bg = "white", device = cairo_pdf)
cat("wrote ", stem, ".{png,pdf}  (n=", nrow(d), ")\n", sep = "")
}
