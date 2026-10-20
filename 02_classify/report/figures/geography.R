# geography.R — unreported co-pathogens in the Pst surveillance compendium. Styling only.
#
# Every number comes from prep_geography.py. Nothing is computed here, including the country
# ordering, which is read from the TSV so the caption and the panel cannot disagree.
#
# TWO PANELS, and they answer different questions, so they are deliberately NOT the same chart:
#   a  how often is each species detected across all 527 libraries. The declared pathogen sits
#      at 100% and is drawn in the declared colour as the positive control; everything below it
#      was unreported. A plain horizontal bar, because the comparison is between species.
#   b  per-country rate for the two most frequent unreported pathogens. PAIRED points joined by
#      a segment rather than dodged bars: the claim is that the two gradients run OPPOSITE ways,
#      and a segment makes the crossing legible where two bar series do not.
#
# Panel b is ordered by Z. tritici rate, which is what makes the opposing gradient visible. The
# order is computed in the prep and carried in the TSV.
#
#   Rscript 02_classify/report/figures/geography.R

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })

HERE <- "02_classify/report/figures"
sp <- read.delim(file.path(HERE, "geography_species.tsv"), stringsAsFactors = FALSE)
gg <- read.delim(file.path(HERE, "geography.tsv"), stringsAsFactors = FALSE,
                 check.names = FALSE)

DECLARED  <- "#8c0404"   # the ring figure's target_cereal red: the declared pathogen
UNREPORT  <- "#22127e"   # host_associate indigo: unreported
ZT        <- "#3A6EA5"   # house wheat blue
PT        <- "#A9512C"   # house sorghum rust-brown, reused here for the second rust

base <- theme_classic(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.8),
        axis.ticks       = element_line(colour = "black", linewidth = 0.8),
        axis.ticks.length = unit(3.5, "pt"),
        axis.text        = element_text(colour = "black"),
        plot.title       = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle    = element_text(size = 10, colour = "grey30"),
        legend.position  = "none",
        plot.margin      = margin(10, 14, 8, 10))

# ── panel a ─────────────────────────────────────────────────────────────────────────────────
# Pathogen guild only. The resident community (Cladosporium allicinum, Dioszegia, Sporobolomyces)
# is detected at comparable coverage and is discussed in the caption, but mixing residents into
# a panel about unreported PATHOGENS would invite the reading that they are co-infections.
# grepl is CASE SENSITIVE and FUNGuild writes "Plant Pathogen". Matching "pathogen" returned
# only the declared species and silently emptied the panel, which looked like a data problem
# rather than a regex one.
a <- sp[sp$n_runs >= 30 &
          (sp$is_declared == "yes" | grepl("Plant Pathogen", sp$guild, fixed = TRUE)), ]
a <- a[order(a$pct_runs), ]
a$species <- factor(a$species, levels = a$species)
a$kind <- ifelse(a$is_declared == "yes", "declared", "unreported")

pA <- ggplot(a, aes(species, pct_runs, fill = kind)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.8) +
  geom_text(aes(label = percent(pct_runs, accuracy = 0.1)),
            hjust = -0.15, size = 3.9, colour = "black") +
  scale_fill_manual(values = c(declared = DECLARED, unreported = UNREPORT)) +
  scale_y_continuous(labels = percent, limits = c(0, 1.18),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  coord_flip() +
  labs(x = NULL, y = "libraries with the species detected",
       title = "a   detected across all 527 libraries",
       subtitle = "red = the pathogen the study reported; indigo = unreported") +
  base + theme(axis.text.y = element_text(face = "italic"))

# ── panel b ─────────────────────────────────────────────────────────────────────────────────
b <- gg[gg$drawn == "yes", c("country", "n_runs",
                             "Zymoseptoria_tritici", "Puccinia_triticina")]
b <- b[order(b$Zymoseptoria_tritici), ]
b$country <- factor(b$country, levels = b$country)
b$lab <- sprintf("%s (%d)", as.character(b$country), b$n_runs)
b$lab <- factor(b$lab, levels = b$lab)

pB <- ggplot(b) +
  geom_segment(aes(x = lab, xend = lab,
                   y = Zymoseptoria_tritici, yend = Puccinia_triticina),
               colour = "grey55", linewidth = 0.9) +
  geom_point(aes(lab, Zymoseptoria_tritici), shape = 21, size = 4.4, stroke = 0.9,
             fill = ZT, colour = "black") +
  geom_point(aes(lab, Puccinia_triticina), shape = 23, size = 4.4, stroke = 0.9,
             fill = PT, colour = "black") +
  scale_y_continuous(labels = percent, limits = c(0, 0.72),
                     breaks = seq(0, 0.6, 0.2), expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(x = NULL, y = "libraries with the species detected",
       title = "b   the two commonest, by country",
       subtitle = "circle = Z. tritici    diamond = P. triticina    (n libraries)") +
  base

p <- pA | pB

ragg::agg_png(file.path(HERE, "geography.png"), width = 4200, height = 1900,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "geography.pdf"), width = 14, height = 6.33)
print(p); invisible(dev.off())
for (nm in c("a", "b")) {
  obj <- if (nm == "a") pA else pB
  ragg::agg_png(file.path(HERE, sprintf("geography_%s.png", nm)),
                width = 1950, height = 1750, res = 300, background = "white")
  print(obj); invisible(dev.off())
}
cat("wrote geography.{png,pdf} and geography_{a,b}.png\n")
