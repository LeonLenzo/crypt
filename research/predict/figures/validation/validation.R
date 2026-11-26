#!/usr/bin/env Rscript
# The prediction-validation figure: four panels, four claims.
#
#   A  fragmentation does not degrade recovery
#   B  predictions are as complete as submitter annotations
#   C  a congener does not cover a species' k-mer space
#   D  large genomes are over-predicted, and that is the limitation
#
# Self-contained mode: axis titles, legends and strip labels kept, since this drops
# straight into the supplement. Run prep.py first.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(scales)
})

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
tsv <- function(f) read.delim(file.path(here, f), stringsAsFactors = FALSE)

# One hue, chosen for the subject: predicted gene sequence. Grey is always the
# comparator or the null, never a third category.
ACCENT <- "#2e7d5b"
LIGHT  <- "#a8d5c2"
GREY   <- "#8d8d8d"

house <- theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey88"),
    axis.line        = element_line(linewidth = 0.8, colour = "grey20"),
    axis.ticks       = element_line(linewidth = 0.8, colour = "grey20"),
    strip.text       = element_text(size = 15),
    plot.tag         = element_text(size = 20, face = "bold"),
    legend.position  = "none",
    plot.margin      = margin(8, 16, 4, 10),
    axis.title.x     = element_text(margin = margin(t = 6)),
    axis.title.y     = element_text(margin = margin(r = 6)))

ital <- function(x) parse(text = sprintf('italic("%s")', x))
abbr <- function(x) ifelse(grepl("^\\w+ sp\\.?$", x), x, sub("^(\\w)\\w+ ", "\\1. ", x))

# ── A ─────────────────────────────────────────────────────────────────────────
a <- tsv("panel_a_benchmark.tsv")
lv <- c("Complete Genome", "Chromosome", "Scaffold", "Contig")
a$assembly_level <- factor(a$assembly_level, levels = lv, labels = c("Complete", lv[-1]))
a$metric <- factor(a$metric, levels = c("recall", "precision"), labels = c("Recall", "Precision"))

pA <- ggplot(a, aes(assembly_level, value)) +
  geom_boxplot(fill = "lightgray", colour = "black", alpha = 0.3,
               linewidth = 1, width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(shape = metric), fill = ACCENT, colour = "black",
              size = 3, stroke = 0.9, alpha = 0.7, width = 0.22, height = 0) +
  facet_wrap(~metric) +
  scale_shape_manual(values = c(Recall = 21, Precision = 22), guide = "none") +
  scale_y_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "Assembly level of the target", y = "k-mer agreement with true CDS") +
  house + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ── B ─────────────────────────────────────────────────────────────────────────
b <- tsv("panel_b_production.tsv")
ref <- tsv("panel_b_reference.tsv")
b <- b[!is.na(b$size_dev_log2), ]
b$flag <- b$busco_complete < 90
lab <- b[b$flag, ]
lab$sp  <- abbr(lab$taxon)
lab$det <- sprintf("%s contigs", comma(lab$contig_count))
lab$hj  <- ifelse(lab$size_dev_log2 > 0.4, 1.12, -0.12)   # keep every label inside the panel
lab$nx  <- ifelse(lab$hj < 0, 1, -1) * sqrt(lab$contig_count) / 900  # clear the marker itself
lab$ny  <- ifelse(lab$contig_count > 10000, -5, 0)        # the 23,880-contig marker is huge

pB <- ggplot(b, aes(size_dev_log2, busco_complete)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ref$q1, ymax = ref$q3,
           fill = GREY, alpha = 0.18) +
  geom_hline(yintercept = ref$median, linewidth = 1, colour = GREY, linetype = "22") +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey75") +
  geom_point(aes(size = contig_count), fill = ACCENT, colour = "black",
             shape = 21, stroke = 0.9, alpha = 0.7) +
  geom_text(data = lab, aes(label = sp, hjust = hj, x = size_dev_log2 + nx,
                            y = busco_complete + ny), vjust = -0.15,
            size = 3.5, fontface = "italic", colour = "grey25") +
  geom_text(data = lab, aes(label = det, hjust = hj, x = size_dev_log2 + nx,
                            y = busco_complete + ny), vjust = 1.35,
            size = 3.3, colour = "grey45") +
  scale_size_continuous(range = c(2, 9), trans = "sqrt", name = "Contigs",
                        breaks = c(100, 2000, 24000), labels = comma) +
  scale_x_continuous(limits = c(-1.75, 1.3), breaks = log2(c(0.5, 1, 2)),
                     labels = c("0.5x", "1x", "2x")) +
  scale_y_continuous(limits = c(20, 101), breaks = seq(20, 100, 20)) +
  labs(x = "Assembly size relative to the genus median",
       y = "BUSCO completeness (%)") +
  house + theme(legend.position = c(0.97, 0.30), legend.justification = c(1, 0.5),
                legend.direction = "horizontal", legend.background = element_blank(),
                legend.title = element_text(size = 13))

# ── C ─────────────────────────────────────────────────────────────────────────
cc <- tsv("panel_c_novelty.tsv")
cc$grp <- factor(ifelse(cc$genus_in_dbv4 == "yes", "yes", "no"), levels = c("yes", "no"),
                 labels = c("Genus already\nin the database", "Genus absent\nfrom the database"))

pC <- ggplot(cc, aes(grp, novelty)) +
  geom_boxplot(fill = "lightgray", colour = "black", alpha = 0.3,
               linewidth = 1, width = 0.55, outlier.shape = NA) +
  geom_jitter(fill = ACCENT, colour = "black", shape = 21,
              size = 3.4, stroke = 0.9, alpha = 0.7, width = 0.2, height = 0) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Predicted k-mers absent from db_v4") +
  house

# ── D ─────────────────────────────────────────────────────────────────────────
d <- tsv("panel_d_overprediction.tsv")
dl <- d[d$overprediction > 2, ]
dl$sp  <- abbr(dl$taxon)
dl$det <- sprintf("precision %.2f", dl$precision)

pD <- ggplot(d, aes(genome_size_mb, overprediction)) +
  geom_hline(yintercept = 1, linewidth = 0.5, colour = "grey75") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, fullrange = FALSE,
              colour = GREY, linewidth = 1, linetype = "22") +
  geom_point(aes(fill = precision), colour = "black", shape = 21,
             size = 4, stroke = 0.9, alpha = 0.85) +
  geom_text(data = dl, aes(label = sp), hjust = 1.12, vjust = -0.75,
            size = 3.6, fontface = "italic", colour = "grey25") +
  geom_text(data = dl, aes(label = det), hjust = 1.12, vjust = 0.75,
            size = 3.4, colour = "grey45") +
  scale_fill_gradient(low = "white", high = ACCENT, limits = c(0.4, 1),
                      breaks = c(0.5, 0.75, 1), name = "Precision") +
  scale_x_continuous(breaks = seq(0, 150, 50), limits = c(0, 155),
                     labels = function(x) paste0(x, " Mb")) +
  scale_y_continuous(limits = c(0.7, 4.15), breaks = 1:4,
                     labels = function(y) paste0(y, "x")) +
  labs(x = "Genome size", y = "Predicted / true CDS bases") +
  house + theme(legend.position = c(0.04, 0.88), legend.justification = c(0, 1),
                legend.direction = "horizontal", legend.background = element_blank(),
                legend.title = element_text(size = 13),
                legend.key.width = unit(1.8, "lines"), legend.key.height = unit(0.7, "lines"))

fig <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "a", theme = house)

agg_png(file.path(here, "validation.png"), width = 3300, height = 2600,
        res = 230, background = "white")
print(fig)
invisible(dev.off())
cat("wrote validation.png\n")
