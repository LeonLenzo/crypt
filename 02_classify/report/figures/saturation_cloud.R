#!/usr/bin/env Rscript
# saturation_cloud.R — one panel per SAMPLE: which of its calls behave like a real organism.
#
# The 03_kraken charts did this per SPECIES, faceted by host, and read the signal off an
# accumulation line. That works there because every point in a panel is the same taxon and so
# shares one ceiling. Per sample it cannot, because every point is a different taxon.
#
# Normalising both axes by the taxon's own db_v4 total makes the null UNIVERSAL:
#     x = reads * kmers_per_spot / M     redundancy available
#     y = distinct / M                   coverage observed
#     null: y = 1 - exp(-x)              one curve, every taxon, any genome size
# so a 2.9 Mb Colletotrichum and a 33 Mb wheat sit on the same scale and the same curve.
#
#   Rscript 02_classify/report/figures/saturation_cloud.R

suppressPackageStartupMessages({ library(ggplot2); library(scales) })

HERE <- "02_classify/report/figures"
d <- read.delim(file.path(HERE, "saturation_cloud.tsv"), stringsAsFactors = FALSE)
d <- d[is.finite(d$redundancy) & is.finite(d$coverage) & d$coverage > 0, ]

COV_FLOOR <- 0.01

# SATURATION CANNOT DISCRIMINATE AT LOW REDUNDANCY, and the figure has to say so. For small x,
# 1 - exp(-x) ~= x, so coverage ~= redundancy MECHANICALLY: every k-mer a read carries is new,
# whether the taxon is real or a handful of cross-mapped loci. Measured on this data the median
# saturation falls 0.67 -> 0.23 -> 0.053 -> 0.019 -> 0.002 across redundancy decades, which is
# the mechanism and not biology. 70% of taxa sit below 1e-3 and are unjudgeable.
#
# That is not a flaw so much as a boundary: the taxa needing adjudication are precisely the
# high-read ones (Oryza at 99,509 reads), and a taxon with too few reads to be judged also has
# too few to clear the coverage floor. The shaded band marks where the statistic is silent.
DISCRIM_X <- 1e-3

# The null, drawn over the full x range rather than fitted. exp() underflows above ~700, which
# is well inside the range a dominant organism reaches, so x is clamped for the curve only.
xs <- 10^seq(log10(min(d$redundancy[d$redundancy > 0])) - 0.3,
             log10(max(d$redundancy)) + 0.3, length.out = 400)
null <- data.frame(x = xs, y = 1 - exp(-pmin(700, xs)))

d$state <- ifelse(d$clears == "yes", "clears both floors",
                  ifelse(d$reads >= 1000, "many reads, fallen off", "trace"))
d$state <- factor(d$state, levels = c("clears both floors", "many reads, fallen off", "trace"))

# Label only what a reader needs to find: the calls that clear, plus the worst offenders.
lab <- do.call(rbind, lapply(split(d, d$run), function(x) {
  keep <- x[x$clears == "yes", ]
  fallen <- x[x$state == "many reads, fallen off", ]
  if (nrow(fallen)) keep <- rbind(keep, fallen[order(-fallen$reads), ][1:min(3, nrow(fallen)), ])
  keep
}))

p <- ggplot(d, aes(redundancy, coverage)) +
  annotate("rect", xmin = min(xs), xmax = DISCRIM_X, ymin = 0, ymax = Inf,
           fill = "grey88", alpha = 0.55) +
  annotate("text", x = sqrt(min(xs) * DISCRIM_X), y = 0.45, label = "saturation\nuninformative",
           size = 2.5, colour = "grey35", lineheight = 0.95) +
  geom_line(data = null, aes(x, y), inherit.aes = FALSE,
            colour = "grey35", linewidth = 0.8) +
  geom_hline(yintercept = COV_FLOOR, linetype = "dashed",
             linewidth = 0.6, colour = "#8c0404") +
  geom_point(aes(fill = state, size = state), shape = 21, stroke = 0.4,
             colour = "black", alpha = 0.85) +
  ggrepel::geom_text_repel(data = lab, aes(label = species), size = 2.3, seed = 1,
                           max.overlaps = 14, min.segment.length = 0,
                           segment.colour = "grey60", segment.size = 0.3) +
  facet_wrap(~ run, ncol = 4) +
  scale_fill_manual(values = c("clears both floors" = "#22127e",
                               "many reads, fallen off" = "#8c0404",
                               "trace" = "grey75"), name = NULL) +
  scale_size_manual(values = c("clears both floors" = 2.6,
                               "many reads, fallen off" = 2.2,
                               "trace" = 1.3), guide = "none") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  # Explicit decade breaks: percent_format(accuracy = 0.001) rendered every label with three
  # decimals, so the dominant organism read "100.000%".
  scale_y_log10(breaks = 10^(-5:0),
                labels = c("0.001%", "0.01%", "0.1%", "1%", "10%", "100%")) +
  annotation_logticks(sides = "bl", size = 0.3,
                      short = unit(2, "pt"), mid = unit(3, "pt"), long = unit(4, "pt")) +
  labs(x = "Redundancy available  (reads × k-mers per spot / taxon's db_v4 minimizers)",
       y = "Coverage observed  (distinct / taxon's db_v4 minimizers)") +
  theme_classic(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.7),
        axis.ticks       = element_line(colour = "black", linewidth = 0.7),
        axis.text        = element_text(colour = "black"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text       = element_text(face = "bold", size = 10),
        legend.position  = "top",
        plot.margin      = margin(10, 14, 8, 10))

n <- length(unique(d$run))
ncol <- min(4, n); nrow <- ceiling(n / ncol)
ragg::agg_png(file.path(HERE, "saturation_cloud.png"),
              width = 1150 * ncol, height = 1050 * nrow + 180, res = 230, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "saturation_cloud.pdf"),
          width = 5 * ncol, height = 4.6 * nrow + 0.8)
print(p); invisible(dev.off())
cat("wrote saturation_cloud.{png,pdf} for", n, "runs\n")
