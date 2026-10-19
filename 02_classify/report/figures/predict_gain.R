#!/usr/bin/env Rscript
# predict_gain.R — what gene prediction would add to db_v4, by route. For a talk.
#
# BAR HEIGHT IS CDS ANNOTATIONS, not species: assemblies carrying a gene model, actual against
# potential. Species count was the first version and Leon's call was that it is not
# informative enough — the database is built from annotated ASSEMBLIES, so that is the honest
# unit for "what would we gain". Species counts stay in predict_gain.tsv and the caption.
#
# Companion to the db_contents ring and deliberately its opposite in form: the ring is a log
# scale over 1,544 tips and answers "what is in db_v4"; this is a plain linear bar chart over
# five categories and answers "what is missing and could be reached".
#
# LINEAR, NOT LOG, and the bars are DODGED rather than stacked.
#   * linear because a bar's LENGTH is its value, which is the one thing a bar chart promises.
#     On a log axis a bar twice as long is more than ten times bigger, and an audience reading
#     it in a talk will not make that adjustment.
#   * the only thing pushing toward log was the `host` route at 6 species, and host has ZERO
#     addable, so it carries no part of the argument. prep_predict_gain.py drops it from the
#     chart and keeps it in the table, which brings the range to 70-1,076 and makes linear
#     comfortable.
#   * dodged because the comparison is current-versus-potential, not a total. Stacking would
#     put the eye on the sum, which is not the claim.
#
# Colours are the ring's own ROUTE_COLOUR, kept in sync by hand here but assigned by the same
# imported precedence rule in the prep, so a species is the same colour in both figures.
#
#   Rscript 02_classify/report/figures/predict_gain.R

suppressPackageStartupMessages({ library(ggplot2); library(scales) })

HERE <- "02_classify/report/figures"
d <- read.delim(file.path(HERE, "predict_gain.tsv"), stringsAsFactors = FALSE)
d <- d[d$drawn == "yes", ]

# Same hexes as db_contents.R.
ROUTE_COLOUR <- c(
  target_cereal   = "#8c0404",
  target_other    = "#e77e7e",
  host_associate  = "#22127e",
  plant_associate = "#add0e4",
  congener        = "#dfd8d8")
ROUTE_WRAP <- c(
  target_cereal   = "pathogen on\nwheat, rice,\nmaize, sorghum",
  target_other    = "plant pathogen,\nnot on those\ncrops",
  host_associate  = "host associate:\nphytobiome\nstudy",
  plant_associate = "plant associate:\nrecorded on\na plant",
  congener        = "congener:\nshares a genus\nwith a pathogen")

ord <- names(ROUTE_COLOUR)
d$route <- factor(d$route, levels = ord)
d$status <- factor(d$status, levels = c("in_db", "addable"),
                   labels = c("in db_v4", "reachable by prediction"))

# Pale = the same hue mixed toward white, so the pair reads as one route at two states rather
# than as two unrelated colours.
lighten <- function(hex, f = 0.62) {
  m <- grDevices::col2rgb(hex)
  grDevices::rgb(t(m + (255 - m) * f), maxColorValue = 255)
}
# Two extra fills carry the LEGEND only. Fill encodes route-and-status, which is ten keys —
# far too many to show — but the convention a viewer needs is just solid-versus-pale. So two
# neutral greys are added as unused factor levels and are the only entries the guide draws.
fills <- c(ROUTE_COLOUR, setNames(lighten(ROUTE_COLOUR), paste0(ord, "_pale")),
           legend_solid = "#4d4d4d", legend_pale = lighten("#4d4d4d"))
d$fillkey <- factor(ifelse(d$status == "in db_v4", as.character(d$route),
                           paste0(as.character(d$route), "_pale")),
                    levels = c(ord, paste0(ord, "_pale"), "legend_solid", "legend_pale"))

p <- ggplot(d, aes(route, cds_annotations, fill = fillkey, group = status)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           colour = "black", linewidth = 1) +
  geom_text(aes(label = comma(cds_annotations)), position = position_dodge(width = 0.78),
            vjust = -0.45, size = 5.2, colour = "black") +
  scale_fill_manual(values = fills, guide = "none") +
  # KEY DRAWN BY HAND. fill encodes route-and-status, which is ten keys, and the usual tricks
  # for showing only two of them both failed: scale breaks on unused factor levels draw the
  # labels with EMPTY glyphs, and guide_legend(override.aes = ...) did not paint them either
  # under ggplot2 4.0.3. Two annotate() rects are exact, need no guide machinery, and put the
  # key where it belongs on a slide.
  annotate("rect", xmin = 3.42, xmax = 3.60, ymin = 1640, ymax = 1730,
           fill = "#4d4d4d", colour = "black", linewidth = 1) +
  annotate("text", x = 3.68, y = 1685, label = "in db_v4 now",
           hjust = 0, size = 5.4, colour = "black") +
  annotate("rect", xmin = 3.42, xmax = 3.60, ymin = 1510, ymax = 1600,
           fill = lighten("#4d4d4d"), colour = "black", linewidth = 1) +
  annotate("text", x = 3.68, y = 1555, label = "reachable by gene prediction",
           hjust = 0, size = 5.4, colour = "black") +
  scale_x_discrete(labels = ROUTE_WRAP) +
  scale_y_continuous(labels = comma, breaks = seq(0, 1600, 200),
                     expand = expansion(mult = c(0, 0.10))) +
  coord_cartesian(ylim = c(0, 1750)) +
  labs(x = NULL, y = "CDS annotations") +
  theme_classic(base_size = 17) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 1),
        axis.ticks       = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(4, "pt"),
        axis.text        = element_text(colour = "black"),
        axis.text.x      = element_text(size = 13, lineheight = 0.95),
        legend.position  = "none",
        plot.margin      = margin(10, 14, 8, 8))

ragg::agg_png(file.path(HERE, "predict_gain.png"), width = 3000, height = 1900,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "predict_gain.pdf"), width = 10, height = 6.33)
print(p); invisible(dev.off())
cat("wrote predict_gain.{png,pdf}\n")
