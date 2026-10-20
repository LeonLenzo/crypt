#!/usr/bin/env Rscript
# detect.R — what 05_detect decided. Decision surface, per-sample counts, species called.
#
# Panel a is the DECISION SURFACE itself and not a picture of one: the axes are the two
# quantities the criterion is computed from, and the two drawn boundaries ARE the two floors.
# Shading is the fraction of each bin that was called, so the boundary reads as a gradient
# instead of collapsing to whichever call happens to dominate a boundary bin.
#
#   Rscript 02_classify/report/figures/detect.R

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })

HERE <- "02_classify/report/figures"
sur <- read.delim(file.path(HERE, "detect_surface.tsv"), stringsAsFactors = FALSE)
ps  <- read.delim(file.path(HERE, "detect_persample.tsv"), stringsAsFactors = FALSE)
spp <- read.delim(file.path(HERE, "detect_species.tsv"), stringsAsFactors = FALSE)

COV_FLOOR <- 0.01
SAT_FLOOR <- 0.05
CROP <- c("Wheat" = "#3A6EA5", "Maize" = "#D9A441", "Rice" = "#45896A",
          "Sorghum" = "#A9512C", "Barley" = "#8C6D8F", "Other cereal" = "#8D9EA9")
GUILD <- c("plant pathogen" = "#8c0404", "endophyte" = "#22127e",
           "saprotroph" = "#45896A", "animal pathogen" = "#8D9EA9",
           "other" = "#D9A441", "unlabelled" = "grey75")

base <- theme_classic(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.8),
        axis.ticks       = element_line(colour = "black", linewidth = 0.8),
        axis.text        = element_text(colour = "black"),
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 10, colour = "grey30"),
        plot.margin      = margin(10, 14, 8, 10))

# ── a. the decision surface ─────────────────────────────────────────────────────────────────
# The shape boundary is coverage = SAT_FLOOR * (1 - exp(-redundancy)), which on log-log axes
# is a curve that flattens to the saturation floor at high redundancy and runs parallel to the
# diagonal at low redundancy -- where, as the rarefaction work showed, it cannot discriminate.
xs <- 10^seq(min(sur$x) - 0.2, max(sur$x) + 0.2, length.out = 400)
bound <- data.frame(x = log10(xs),
                    y = log10(SAT_FLOOR * (1 - exp(-pmin(700, xs)))))
bound <- bound[is.finite(bound$y), ]

pa <- ggplot(sur, aes(x, y)) +
  geom_tile(aes(fill = frac_called, alpha = log10(n)), width = 0.062, height = 0.062) +
  geom_hline(yintercept = log10(COV_FLOOR), linetype = "dashed",
             linewidth = 0.9, colour = "black") +
  geom_line(data = bound, aes(x, y), inherit.aes = FALSE,
            linewidth = 0.9, colour = "black") +
  annotate("text", x = min(sur$x) + 0.3, y = log10(COV_FLOOR) + 0.28,
           label = "abundance floor: coverage 1%", hjust = 0, size = 3, colour = "black") +
  annotate("text", x = max(sur$x) - 0.2, y = log10(SAT_FLOOR) - 0.75,
           label = "shape test:\nsaturation 0.05", hjust = 1, size = 3,
           colour = "black", lineheight = 0.95) +
  scale_fill_gradientn(colours = c("grey88", "#f4c9c9", "#c94c4c", "#8c0404"),
                       limits = c(0, 1), labels = percent_format(1),
                       name = "fraction\ncalled") +
  scale_alpha_continuous(range = c(0.35, 1), guide = "none") +
  scale_x_continuous(labels = function(v) parse(text = paste0("10^", v))) +
  scale_y_continuous(labels = function(v) parse(text = paste0("10^", v))) +
  labs(title = "a  The decision surface",
       subtitle = "300,061 candidates; shading is the fraction called",
       x = "Redundancy available", y = "Coverage observed") +
  base + theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.height = unit(10, "pt"), legend.key.width = unit(8, "pt"),
               legend.title = element_text(size = 8), legend.text = element_text(size = 8))

# ── b. pathogens per sample ─────────────────────────────────────────────────────────────────
ps$crop <- factor(ps$crop, levels = names(CROP))
ps <- ps[!is.na(ps$crop), ]
ps$np <- pmin(ps$n_pathogens, 8)
agg <- as.data.frame(table(crop = ps$crop, n = ps$np))
agg <- agg[agg$Freq > 0, ]
pb <- ggplot(agg, aes(as.integer(as.character(n)), Freq, fill = crop)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.82) +
  geom_vline(xintercept = 1.5, linetype = "dotted", linewidth = 0.7, colour = "grey30") +
  annotate("text", x = 1.62, y = Inf, vjust = 1.5, hjust = 0, size = 3, colour = "grey25",
           label = "co-infected →") +
  scale_fill_manual(values = CROP, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = 0:8, labels = c(0:7, "8+")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.06))) +
  labs(title = "b  Plant pathogens per sample",
       subtitle = "guild labelled after detection, never filtered before",
       x = "Plant pathogens called", y = "Runs") +
  base + theme(legend.position = c(0.98, 0.72), legend.justification = c(1, 1),
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.size = unit(10, "pt"), legend.text = element_text(size = 8.5))

# ── c. what is being called ─────────────────────────────────────────────────────────────────
top <- head(spp[order(-spp$n_runs), ], 22)
top$species <- factor(top$species, levels = rev(top$species))
top$guild <- factor(top$guild, levels = names(GUILD))
pc <- ggplot(top, aes(n_runs, species, fill = guild)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.74) +
  scale_fill_manual(values = GUILD, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "c  Which species are called",
       subtitle = "runs carrying each, strongest guild claim shown",
       x = "Runs", y = NULL) +
  base + theme(axis.text.y = element_text(face = "italic", size = 8.5),
               legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
               legend.background = element_rect(fill = "white", colour = NA),
               legend.key.size = unit(10, "pt"), legend.text = element_text(size = 8.5))

p <- pa | pb | pc

ragg::agg_png(file.path(HERE, "detect.png"), width = 4900, height = 1600,
              res = 240, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "detect.pdf"), width = 17, height = 5.6)
print(p); invisible(dev.off())
for (nm in c("a", "b", "c")) {
  ragg::agg_png(file.path(HERE, sprintf("detect_%s.png", nm)),
                width = 1700, height = 1600, res = 240, background = "white")
  print(switch(nm, a = pa, b = pb, c = pc)); invisible(dev.off())
}
cat("wrote detect.{png,pdf} and detect_{a,b,c}.png\n")
