#!/usr/bin/env Rscript
# rarefaction.R — what library depth does to a detection. Three panels, left to right:
# mechanism, consequence, cost.
#
# LABELS ARE ON. This figure is for working out what the result means, not for assembly in
# Illustrator, so panel titles, axis titles and the project annotations all stay. The other
# 02_classify figures strip them because their panels get placed by hand; say the word and
# there is an --assembly mode to add.
#
# EVERY AXIS IS Gkmer, NOT Gb. At k=35 a read of length L yields L-34 k-mers, so k-mer yield
# per base is (L-34)/L: 0.33 at 51 nt, 0.77 at 150 nt. Read length is worth 2.3x independent
# of depth, so Gb alone would put PRJDB7234 and a 150 nt library of the same size in the same
# place when they are a factor of 2.3 apart. The measured result is that Gkmer is SUFFICIENT:
# at matched Gkmer the 51 nt arm retains as many detections as the native one.
#
#   Rscript 02_classify/report/figures/rarefaction.R

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

HERE <- "02_classify/report/figures"
cur <- read.delim(file.path(HERE, "rarefaction_curve_plot.tsv"), stringsAsFactors = FALSE)
ret <- read.delim(file.path(HERE, "rarefaction_retained.tsv"), stringsAsFactors = FALSE)
cost <- read.delim(file.path(HERE, "rarefaction_threshold_cost.tsv"), stringsAsFactors = FALSE)

COV_FLOOR <- 0.01

# The real projects, so a position on the x axis reads as a statement about the cohort rather
# than about an abstract depth. These two are the whole reason the experiment exists.
# `ypos` staggers the labels: 0.019 and 0.036 are a factor of 1.9 apart, which on a log axis
# spanning 160x is close enough that the two captions collided. `hj` flips the right-hand one
# inward so it is not clipped by the panel edge.
ANCHOR <- data.frame(
  gk    = c(0.019, 0.036, 3.007),
  label = c("PRJDB7234\n1,557 runs, rice", "PRJNA726716\n104 runs", "PRJEB39201\n538 runs"),
  ypos  = c(1.00, 0.34, 1.00),
  hj    = c(-0.07, -0.07, 1.07))

# House crop palette, identical to 01_literature's map figure so a crop is one colour
# across the thesis.
CROP <- c("Wheat" = "#3A6EA5", "Maize" = "#D9A441", "Rice" = "#45896A",
          "Sorghum" = "#A9512C", "Barley" = "#8C6D8F", "Other cereal" = "#8D9EA9")

base <- theme_classic(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.8),
        axis.ticks       = element_line(colour = "black", linewidth = 0.8),
        axis.ticks.length = unit(3.5, "pt"),
        axis.text        = element_text(colour = "black"),
        plot.title       = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle    = element_text(size = 10.5, colour = "grey30"),
        plot.margin      = margin(10, 14, 8, 10))

# ymax positions the project labels at the top of the panel; it has no sensible default
# because each panel's y scale is different (log coverage, a count, a count of runs).
anchor_layer <- function(ymax = NULL, show_labels = TRUE) {
  # inherit.aes = FALSE on both: the enclosing plots map group/colour/shape to columns that
  # exist in the curve data and not in ANCHOR, and a layer silently inheriting them fails at
  # build time with "object 'key' not found" rather than at the call.
  l <- list(geom_vline(data = ANCHOR, aes(xintercept = gk), inherit.aes = FALSE,
                       linetype = "dotted", linewidth = 0.7, colour = "grey45"))
  if (show_labels && !is.null(ymax))
    l <- c(l, list(geom_text(data = ANCHOR, aes(x = gk, y = ymax * ypos, label = label,
                                                hjust = hj),
                             inherit.aes = FALSE,
                             vjust = 1, size = 2.9, lineheight = 0.95, colour = "grey30")))
  l
}

# ── a. mechanism ────────────────────────────────────────────────────────────────────────────
# One line per detection. Coloured by its coverage at the DEEPEST rung, which is the point:
# the abundant detections are saturated and have k-mer reserve to lose, the sparse ones sit on
# the linear part of their accumulation curve and lose k-mer space read for read. So shallow
# sequencing preferentially preserves the organism the study already knew about.
a <- cur[cur$arm == "native", ]
a$key <- paste(a$run, a$species)
deepest <- do.call(rbind, lapply(split(a, a$key), function(d) d[which.max(d$rung_gkmer), ]))
a$full <- deepest$coverage[match(a$key, deepest$key)]
# A coverage of exactly 0 cannot be drawn on a log axis and is not the same as "missing": it
# means no distinct minimizer of that taxon survived. Floored to half the smallest positive
# value so the point is visible AND obviously at the bottom.
floor_val <- min(a$coverage[a$coverage > 0]) / 2
a$cov_plot <- ifelse(a$coverage <= 0, floor_val, a$coverage)

pa <- ggplot(a, aes(rung_gkmer, cov_plot, group = key, colour = full)) +
  anchor_layer(ymax = max(a$cov_plot) * 2.2) +
  geom_hline(yintercept = COV_FLOOR, linetype = "dashed", linewidth = 0.9,
             colour = "#8c0404") +
  annotate("text", x = 0.016, y = COV_FLOOR * 1.35, label = "1% coverage floor",
           hjust = 0, size = 3.1, colour = "#8c0404") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 1.5) +
  scale_colour_viridis_c(option = "viridis", direction = -1, trans = "log10",
                         name = "coverage at\nfull depth", labels = percent_format(1)) +
  scale_x_log10(breaks = c(0.019, 0.1, 0.5, 3), labels = c("0.019", "0.1", "0.5", "3")) +
  # Explicit decade breaks. percent_format() rendered the floored zeros as "0.0%", which reads
  # as a measured zero rather than as "no minimizer of this taxon survived"; the break labels
  # now stop at 0.01% and the floored points sit visibly below the lowest one.
  scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
                labels = c("0.01%", "0.1%", "1%", "10%", "100%")) +
  labs(title = "a  Every detection fades, but not equally",
       subtitle = "22 established detections, 5 wheat libraries",
       x = "Library depth (Gkmer)", y = "Coverage against db_v4") +
  base + theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.height = unit(9, "pt"), legend.key.width = unit(7, "pt"),
               legend.title = element_text(size = 8.5), legend.text = element_text(size = 8))

# ── b. consequence ──────────────────────────────────────────────────────────────────────────
# The count is the point: co-infection COUNT is the dependent variable, so the quantity that
# matters is how many of a sample's pathogens remain visible, not how any one of them scores.
# `shallow` is drawn to show it lands on top of `native`: read length costs nothing beyond
# Gkmer, which was not the expected result and is why one threshold on one number is enough.
ret$arm_lab <- factor(ret$arm, levels = c("native", "shallow"),
                      labels = c("as sequenced", "51 nt single-end"))
# FRACTION OF THOSE MEASURED, not a raw count, and this is a correctness fix rather than a
# style choice. Not every run reaches every rung -- a run can only be rarefied DOWN, so the
# 3 Gkmer rung is carried by 3 of the 5 runs and the shallow arm's 0.5 rung by 10 of 22
# detections. Plotting the count made retention appear to FALL at the deep end, when what fell
# was the denominator. As a fraction both arms are monotone, which is what the data says.
ret$frac <- ret$retained_both / ret$detections_measured
ret$n_lab <- ret$detections_measured
pb <- ggplot(ret, aes(rung_gkmer, frac, colour = arm_lab, shape = arm_lab)) +
  anchor_layer(show_labels = FALSE) +
  geom_hline(yintercept = 18 / 22, linetype = "dashed", linewidth = 0.8, colour = "grey40") +
  annotate("text", x = 0.0165, y = 18 / 22 + 0.035, hjust = 0, size = 3.1, colour = "grey30",
           label = "ceiling: 18 of 22 recoverable at any depth") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.9, fill = "white", stroke = 0.9) +
  geom_text(aes(label = n_lab), vjust = 2.1, size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c("as sequenced" = "#22127e",
                                 "51 nt single-end" = "#D9A441"), name = NULL) +
  scale_shape_manual(values = c("as sequenced" = 21, "51 nt single-end" = 24), name = NULL) +
  scale_x_log10(breaks = c(0.019, 0.1, 0.5, 3), labels = c("0.019", "0.1", "0.5", "3")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = percent_format(1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "b  How many of a sample's pathogens stay visible",
       # Short enough not to run into panel c's title at this width. What the numerals mean
       # is in the caption, which is where a reader looks once rather than every time.
       subtitle = "numerals: detections measurable at that rung",
       x = "Library depth (Gkmer)", y = "Detections clearing both floors") +
  base + theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.height = unit(12, "pt"), legend.text = element_text(size = 8.5))

# ── c. cost ─────────────────────────────────────────────────────────────────────────────────
# Stacked, because the question is "how many runs do I lose in total and who pays", and the
# answer is that rice pays almost all of it until 0.5 Gkmer, where maize starts to.
cost$crop <- factor(cost$crop, levels = names(CROP))
# Candidate cut-offs only. 2 and 3 Gkmer are rungs of the EXPERIMENT, not thresholds anyone
# would adopt -- at 2 Gkmer the cohort loses 4,000 runs including most of wheat -- and leaving
# them in let two implausible bars set the y scale and flatten the range that matters.
cost <- cost[!is.na(cost$crop) & cost$lost > 0 & cost$rung_gkmer <= 1.0, ]
pc <- ggplot(cost, aes(factor(rung_gkmer), lost, fill = crop)) +
  geom_col(colour = "black", linewidth = 0.4, width = 0.78) +
  scale_fill_manual(values = CROP, name = NULL, drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.06))) +
  labs(title = "c  A depth threshold is a rice filter",
       subtitle = "cohort runs dropped at each cut-off",
       x = "Threshold (Gkmer)", y = "Runs dropped") +
  base + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.size = unit(10, "pt"), legend.text = element_text(size = 8.5))

p <- pa | pb | pc

ragg::agg_png(file.path(HERE, "rarefaction.png"), width = 4800, height = 1600,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "rarefaction.pdf"), width = 16, height = 5.33)
print(p); invisible(dev.off())

for (nm in c("a", "b", "c")) {
  obj <- switch(nm, a = pa, b = pb, c = pc)
  ragg::agg_png(file.path(HERE, sprintf("rarefaction_%s.png", nm)),
                width = 1700, height = 1600, res = 300, background = "white")
  print(obj); invisible(dev.off())
}
cat("wrote rarefaction.{png,pdf} and rarefaction_{a,b,c}.png\n")
