#!/usr/bin/env Rscript
# plot_resolution.R — draw the purpose test from what prep_resolution.py decided.
#
#   python release/03_classify/figures/resolution/prep_resolution.py   # decides
#   Rscript release/03_classify/figures/resolution/plot_resolution.R   # draws
#
# Reads  resolution_curves.tsv   one point per run per bucket, already transformed
#        resolution_ticks.tsv    where the axis labels go, in the same space
#        resolution_summary.tsv  the medians and the counts, for both panels
# Writes resolution_reads.png    the headline in three numbers
#        resolution.png          the per-run evidence behind it
#
# THIS SCRIPT DECIDES NOTHING. The symmetric-log transform is applied in prep, so the x
# values and the tick positions arrive in one space and R cannot place a label by a
# different rule than it placed the data.
#
# ## What the reader is meant to do with it
#
# Find the vertical zero line and see which side of it each curve lies on. A curve
# entirely to the right of zero means every run gained reads at that rank. The median is
# where a curve crosses 0.5, and the distance between the species and genus medians is the
# size of the shift in the typical run; the tails are the few runs where it is large.
#
# ## House style
#
# leon-figure-style: clean lines, little text, white background set explicitly. The curves
# are the only ink that matters, so the panel has one axis line each way and no box.

suppressPackageStartupMessages({library(ggplot2); library(grid)})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(HERE) || HERE == "") HERE <- "release/03_classify/figures/resolution"
rd <- function(f) read.delim(file.path(HERE, f), sep = "\t", quote = "",
                             check.names = FALSE, colClasses = "character")

cur <- rd("resolution_curves.tsv")
tic <- rd("resolution_ticks.tsv")
sm  <- rd("resolution_summary.tsv")
cur$x <- as.numeric(cur$x); cur$y <- as.numeric(cur$y)
tic$x <- as.numeric(tic$x)
sm$median_delta_pp <- as.numeric(sm$median_delta_pp)

# Buckets keep prep's order, which is the rank ladder: species, genus, above genus.
cur$bucket <- factor(cur$bucket, levels = sm$bucket)
pal <- setNames(sm$colour, sm$bucket)

# The curve is a STEP, not a line: an ECDF over 1,560 runs is a step function and drawing
# it as a smooth line would imply values between the runs.
p <- ggplot(cur, aes(x = x, y = y, colour = bucket)) +
  geom_vline(xintercept = 0, colour = "#4a4a4a", linewidth = 0.5, linetype = "22") +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(values = pal, name = NULL,
                      guide = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2))) +
  scale_x_continuous(breaks = tic$x, labels = tic$label,
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "", "0.5", "", "1"),
                     expand = expansion(mult = 0.015)) +
  labs(x = "db2 \u2212 db1p (percentage points)", y = "share of runs") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(colour = "#ececec", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.7),
    axis.ticks       = element_line(colour = "black", linewidth = 0.7),
    axis.ticks.length = unit(3, "pt"),
    axis.title       = element_text(size = 10, colour = "#1a1a1a"),
    axis.text        = element_text(size = 9.5, colour = "#3a3a3a"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9.5, colour = "#1a1a1a",
                                    margin = margin(l = 4, r = 18)),
    legend.margin    = margin(2, 0, 0, 0),
    plot.margin      = margin(8, 14, 6, 6))

SIDE <- 5.0          # both panels square, and the same square
ragg::agg_png(file.path(HERE, "resolution.png"), width = SIDE, height = SIDE,
              units = "in", res = 600, background = "white")
print(p)
invisible(dev.off())

# ── the headline panel: where the reads went, in reads ────────────────────────
# Three numbers, one bar each, because that is the whole result: the predictions take reads
# off species and the genus node receives them. It sits BESIDE the ECDF rather than
# replacing it -- the bar says what happened across the cohort, the curves say that it
# happened in nearly every run and not in a few deep libraries.
#
# THE TWO PANELS COUNT DIFFERENT THINGS, deliberately: reads here, share of a run's fungal
# reads there. A run can gain reads at a rank while that rank's share of its reads falls,
# which is why the counts under these bars are not the counts the curves cross at.
sm$total_delta_reads <- as.numeric(sm$total_delta_reads)
sm$n_up_reads <- as.integer(sm$n_up_reads); sm$n_down_reads <- as.integer(sm$n_down_reads)
sm$n <- as.integer(sm$n)
sm$bucket <- factor(sm$bucket, levels = rev(sm$bucket))     # species at the top
sm$M <- sm$total_delta_reads / 1e6
sm$moved <- ifelse(sm$M >= 0, sm$n_up_reads, sm$n_down_reads)
# The number sits on the far side of the bar from zero, so a label never lies over its bar.
# The run count sits UNDER the bar and is anchored at the zero line rather than at the
# bar's end: anchored at the end it ran off the panel on the longest bar, which is exactly
# the bar a reader most wants the count for.
sm$hj <- ifelse(sm$M >= 0, -0.12, 1.12)
sm$note_hj <- ifelse(sm$M >= 0, 0, 1)
# a TYPOGRAPHIC MINUS, not a hyphen: the hyphen sits high and short beside a numeral
sm$text <- sub("-", "\u2212", sprintf("%+.1f M", sm$M), fixed = TRUE)
sm$note <- sprintf("%s of %s runs", format(sm$moved, big.mark = ","),
                   format(sm$n, big.mark = ","))
LIMIT <- max(abs(sm$M)) * 1.30

q <- ggplot(sm, aes(x = M, y = bucket, fill = colour)) +
  geom_vline(xintercept = 0, colour = "#4a4a4a", linewidth = 0.5) +
  geom_col(width = 0.55) +
  scale_fill_identity() +
  geom_text(aes(label = text, hjust = hj), size = 3.6, colour = "#1a1a1a",
            fontface = "bold") +
  geom_text(aes(x = 0, label = note, hjust = note_hj), size = 3.0, colour = "#707070",
            nudge_y = -0.33, nudge_x = 0) +
  scale_x_continuous(limits = c(-LIMIT, LIMIT), breaks = c(-50, -25, 0, 25, 50),
                     labels = function(x) sub("-", "\u2212", x, fixed = TRUE),
                     expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(add = 0.7)) +
  coord_cartesian(clip = "off") +
  labs(x = "db2 \u2212 db1p (million reads)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.line.x      = element_line(colour = "black", linewidth = 0.7),
    axis.ticks.x     = element_line(colour = "black", linewidth = 0.7),
    axis.ticks.y     = element_blank(),
    axis.ticks.length = unit(3, "pt"),
    axis.title.x     = element_text(size = 10, colour = "#1a1a1a", margin = margin(t = 6)),
    axis.text.x      = element_text(size = 9.5, colour = "#3a3a3a"),
    axis.text.y      = element_text(size = 11, colour = "#1a1a1a", hjust = 1,
                                    margin = margin(r = 6)),
    # the labels are drawn past the end of their bars with clipping off, so the margins
    # hold them rather than the scale: widening the scale instead would shrink every bar
    # to make room for text.
    plot.margin      = margin(10, 44, 6, 6))

ragg::agg_png(file.path(HERE, "resolution_reads.png"), width = SIDE, height = SIDE,
              units = "in", res = 600, background = "white")
print(q)
invisible(dev.off())

cat(sprintf("wrote %s/resolution_reads.png  %.2f x %.2f in\n", HERE, SIDE, SIDE))
cat(sprintf("wrote %s/resolution.png  %.2f x %.2f in\n", HERE, SIDE, SIDE))
for (i in seq_len(nrow(sm)))
  cat(sprintf("  %-12s median %+6.3f pp, up in %s of %s runs\n", sm$bucket[i],
              sm$median_delta_pp[i], format(as.integer(sm$n_up[i]), big.mark = ","),
              format(as.integer(sm$n[i]), big.mark = ",")))
