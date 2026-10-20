#!/usr/bin/env Rscript
# detect_runs.R — the detection decision, one sample per panel.
#
# Same axes as saturation_cloud.R and for the same reason: normalising both by the taxon's own
# db_v4 minimizer total makes y = 1-exp(-x) a single null curve for every taxon, so a 2.9 Mb
# Colletotrichum and 33 Mb wheat are comparable on one panel. What is added here is the CALL.
#
# Points are filled by guild when 05_detect accepted them and left hollow grey when it did not,
# so "which of this sample's organisms are real" is a question about the coloured points only.
#
#   Rscript 02_classify/report/figures/detect_runs.R

suppressPackageStartupMessages({ library(ggplot2); library(scales) })

HERE <- "02_classify/report/figures"
d <- read.delim(file.path(HERE, "detect_runs.tsv"), stringsAsFactors = FALSE)
d <- d[is.finite(d$redundancy) & is.finite(d$coverage) & d$coverage > 0, ]
d$panel <- factor(d$panel, levels = unique(d$panel[order(d$panel_order)]))

COV_FLOOR <- 0.01
SAT_FLOOR <- 0.05
GUILD <- c("plant pathogen" = "#8c0404", "endophyte" = "#22127e",
           "saprotroph" = "#45896A", "animal pathogen" = "#8D9EA9",
           "other" = "#D9A441", "unlabelled" = "#add0e4")

xs <- 10^seq(log10(min(d$redundancy)) - 0.3, log10(max(d$redundancy)) + 0.3,
             length.out = 500)
null  <- data.frame(x = xs, y = 1 - exp(-pmin(700, xs)))
# The shape test as a boundary in these coordinates: a scaled copy of the null.
shape <- data.frame(x = xs, y = SAT_FLOOR * (1 - exp(-pmin(700, xs))))

d$state <- ifelse(d$call == "yes", "called", "rejected")
called <- d[d$call == "yes", ]
# Label the calls, plus anything rejected despite a large read count -- those are the
# decisions a reader is most likely to want to check.
loud <- d[d$call == "no" & d$reads >= 20000, ]
loud <- do.call(rbind, lapply(split(loud, loud$run),
                              function(x) head(x[order(-x$reads), ], 2)))
lab <- rbind(called, loud)

p <- ggplot(d, aes(redundancy, coverage)) +
  geom_line(data = null, aes(x, y), inherit.aes = FALSE,
            colour = "grey40", linewidth = 0.8) +
  geom_line(data = shape, aes(x, y), inherit.aes = FALSE,
            colour = "black", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = COV_FLOOR, linetype = "dashed",
             linewidth = 0.7, colour = "black") +
  geom_point(data = d[d$call == "no", ], shape = 21, size = 1.5, stroke = 0.35,
             colour = "grey55", fill = NA, alpha = 0.6) +
  geom_point(data = called, aes(fill = guild), shape = 21, size = 3.1,
             stroke = 0.6, colour = "black") +
  ggrepel::geom_text_repel(data = lab, aes(label = species, colour = state),
                           size = 2.25, seed = 1, max.overlaps = 11,
                           min.segment.length = 0, segment.size = 0.25,
                           segment.colour = "grey65", show.legend = FALSE) +
  facet_wrap(~ panel, ncol = 4) +
  scale_fill_manual(values = GUILD, name = NULL, drop = FALSE) +
  scale_colour_manual(values = c(called = "black", rejected = "grey45"), guide = "none") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = 10^(-6:0),
                labels = c("0.0001%", "0.001%", "0.01%", "0.1%", "1%", "10%", "100%")) +
  annotation_logticks(sides = "bl", size = 0.3, short = unit(2, "pt"),
                      mid = unit(3, "pt"), long = unit(4, "pt")) +
  labs(x = "Redundancy available  (reads × k-mers per spot / taxon's db_v4 minimizers)",
       y = "Coverage observed  (distinct / taxon's db_v4 minimizers)",
       caption = paste("solid curve: null for reads spread uniformly   ·",
                       "long-dash: shape test (saturation 0.05)   ·",
                       "dash: abundance floor (coverage 1%)")) +
  theme_classic(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.7),
        axis.ticks       = element_line(colour = "black", linewidth = 0.7),
        axis.text        = element_text(colour = "black", size = 8.5),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text       = element_text(face = "bold", size = 9.5),
        legend.position  = "top",
        plot.caption     = element_text(size = 8.5, colour = "grey30", hjust = 0.5),
        plot.margin      = margin(8, 14, 6, 10))

n <- length(unique(d$panel)); ncol <- min(4, n); nrow <- ceiling(n / ncol)
ragg::agg_png(file.path(HERE, "detect_runs.png"),
              width = 1200 * ncol, height = 1080 * nrow + 230, res = 225,
              background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "detect_runs.pdf"),
          width = 5.2 * ncol, height = 4.7 * nrow + 1)
print(p); invisible(dev.off())
cat("wrote detect_runs.{png,pdf} for", n, "runs\n")
