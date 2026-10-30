# depth_funnel.R -- why each crop leaves the analysis, and the sequencing depth behind it.
# Styling only; every number comes from prep_depth_funnel.py.
#
# THE POINT OF THE FIGURE. "Rice is 0% co-infected" would be a result. It is not one. Every
# rice run is excluded before anything is counted, either because the study suppressed disease
# by design or because the library is too shallow to interpret a negative. This draws both
# reasons next to the depth that produces the second.
#
# PANEL A IS ON A LOG AXIS because the crops differ by orders of magnitude, not percentages:
# rice's median is 0.050 Gkmer against wheat's 3.07. On a linear axis rice is a line on the
# left edge and the comparison cannot be read.
#
# WHY Gkmer AND NOT READS. ENA's read_len is bases per SPOT, so a paired library's per-read
# length is half of it, and at k=35 a 51 nt read yields 17 k-mers against 116 for a 150 nt
# read. Rice is 1,636 single-end 51 nt runs of 2,784. Ranking on reads or bases hides that.
#
#   Rscript 02_classify/report/figures/depth_funnel/depth_funnel.R

suppressPackageStartupMessages({
  library(ggplot2); library(scales); library(ragg); library(patchwork)
})

HERE <- "02_classify/report/figures/depth_funnel"
d    <- read.delim(file.path(HERE, "depth_funnel.tsv"), stringsAsFactors = FALSE)
s    <- read.delim(file.path(HERE, "depth_funnel_summary.tsv"), stringsAsFactors = FALSE)

FLOOR <- 0.5

# Chosen for this figure. Retained is the one saturated colour; every exclusion is a reason,
# not a value, so they read as a set without competing with it.
FATE <- c("retained"                                   = "#1b7a5a",
          "below the depth floor"                      = "#b8472f",
          "disease suppressed or manipulated by design" = "#6b5b8a",
          "exposure undetermined"                      = "grey72")
d$fate <- factor(d$fate, levels = names(FATE))

# Crops ordered by median depth, so panel A reads as a gradient and panel B keeps the same rows.
ord <- s$crop[order(s$median_gkmer)]
d$crop <- factor(d$crop, levels = ord)
s$crop <- factor(s$crop, levels = ord)

base <- theme_minimal(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y      = element_text(size = 12, colour = "grey15"),
        plot.title       = element_text(face = "bold", size = 13,
                                        margin = margin(b = 6)),
        legend.position  = "bottom",
        legend.title     = element_blank())

dd <- d[!is.na(d$gkmer) & d$gkmer > 0, ]

a <- ggplot(dd, aes(gkmer, crop)) +
  annotate("rect", xmin = 1e-6, xmax = FLOOR, ymin = -Inf, ymax = Inf,
           fill = "#f8ece8", colour = NA) +
  geom_violin(scale = "width", fill = "grey80", colour = NA, width = 0.85) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "grey15",
               linewidth = 0.6) +
  geom_vline(xintercept = FLOOR, colour = "#b8472f", linewidth = 0.9) +
  # Labels sit in a reserved decade to the right of the data rather than over it. The deepest
  # run in the cohort is 32 Gkmer, so nothing is hidden by extending the axis to 400.
  geom_text(data = s, aes(x = 480, y = crop, label = sprintf("%s%% below", pct_below_floor)),
            inherit.aes = FALSE, hjust = 1, size = 3.9, colour = "#b8472f") +
  # Lower limit is the true non-zero minimum, so nothing is silently clipped. 24 rice runs
  # sit between 1e-5 and 1e-4, which is four orders below the floor and worth seeing.
  scale_x_log10(limits = c(1e-5, 520),
                breaks = c(1e-5, 1e-3, 0.1, 1, 10),
                labels = c("0.00001", "0.001", "0.1", "1", "10")) +
  labs(title = "a   Sequencing depth, as fungal-detectable k-mer space",
       x = "Gkmer per run (log scale)", y = NULL) +
  base

b <- ggplot(d, aes(y = crop, fill = fate)) +
  geom_bar(width = 0.68) +
  scale_fill_manual(values = FATE, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), labels = comma) +
  labs(title = "b   What happens to each run", x = "runs in the curated cohort", y = NULL) +
  base + theme(axis.text.y = element_blank())

out <- (a | b) + plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "Why the rice studies are excluded before anything is counted",
    caption = paste0(
      "The red line is the 0.5 Gkmer inclusion floor; the shaded band is below it. Every one of 2,784 rice runs is excluded: 1,123 because the study sheltered,\n",
      "protected or did not report its disease exposure, 104 because the plants were raised in a chamber before transfer, and the remaining 1,557 on depth.\n",
      "Rice is 1,636 single-end 51 nt runs, which at k=35 yield 17 k-mers per read against 116 for a 150 nt read, so the gap is wider than base counts suggest.\n",
      "A log axis cannot show zero, so panel a omits 45 rice runs with effectively no bases recorded and 2 sorghum runs with no depth metadata. Panel b counts all 7,688."),
    theme = theme(plot.title   = element_text(face = "bold", size = 16,
                                              margin = margin(b = 8)),
                  plot.caption = element_text(size = 9, colour = "grey35", hjust = 0,
                                              margin = margin(t = 8)),
                  plot.background = element_rect(fill = "white", colour = NA),
                  legend.position = "bottom"))

ragg::agg_png(file.path(HERE, "depth_funnel.png"), width = 4200, height = 1900,
              res = 300, background = "white")
print(out); invisible(dev.off())
cairo_pdf(file.path(HERE, "depth_funnel.pdf"), width = 14, height = 6.33)
print(out); invisible(dev.off())
cat("wrote depth_funnel.{png,pdf}\n")
