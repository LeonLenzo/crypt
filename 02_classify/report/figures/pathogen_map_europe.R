# pathogen_map_europe.R — the Europe window of pathogen_map.R. Styling only.
#
# Companion, not a replacement. The world figure answers "is this species global or regional";
# this one answers "where in the region", and it exists because 67% of localities and 69% of
# libraries sit in a window that is about a fifth of each world panel.
#
# TWO DIFFERENCES FROM THE WORLD FIGURE, both earned by the extra space:
#
#   1. EVERY SAMPLED LOCALITY IS DRAWN IN EVERY PANEL, as a filled grey point rather than an
#      open one, and the detections sit on top. At world scale the grey layer is a texture; at
#      this scale a reader can follow an individual site across panels and see that the same
#      field carries septoria and not tan spot. That comparison is the reason to zoom.
#
#   2. EACH PANEL PRINTS ITS EUROPEAN DETECTION RATE, from pathogen_map_europe.tsv. The map
#      shows WHERE and point density follows sampling effort; the printed rate is the quantity
#      with a denominator. Putting it in the panel stops the eye reading density as frequency,
#      which is the single easiest mistake to make with this figure.
#
# The basemap follows the house convention: land grey92, borders in the ocean colour so they
# read as hairline separations rather than competing with the points.
#
#   Rscript 02_classify/report/figures/pathogen_map_europe.R

suppressPackageStartupMessages({
  library(ggplot2); library(maps); library(scales); library(ragg)
})

HERE  <- "02_classify/report/figures"
d     <- read.delim(file.path(HERE, "pathogen_map.tsv"), stringsAsFactors = FALSE)
sites <- read.delim(file.path(HERE, "pathogen_map_sites.tsv"), stringsAsFactors = FALSE)
eur   <- read.delim(file.path(HERE, "pathogen_map_europe.tsv"), stringsAsFactors = FALSE)

DRAW <- c("Zymoseptoria tritici", "Puccinia triticina", "Blumeria graminis",
          "Microdochium nivale", "Parastagonospora nodorum", "Pyrenophora tritici-repentis")
COMMON <- c("Zymoseptoria tritici"         = "septoria tritici blotch",
            "Puccinia triticina"           = "brown rust",
            "Blumeria graminis"            = "powdery mildew",
            "Microdochium nivale"          = "snow mould",
            "Parastagonospora nodorum"     = "septoria nodorum blotch",
            "Pyrenophora tritici-repentis" = "tan spot")
FILL <- c("Zymoseptoria tritici"         = "#3A6EA5",
          "Puccinia triticina"           = "#A9512C",
          "Blumeria graminis"            = "#45896A",
          "Microdochium nivale"          = "#8C6D8F",
          "Parastagonospora nodorum"     = "#D9A441",
          "Pyrenophora tritici-repentis" = "#8c0404")

XLIM <- c(-11, 32); YLIM <- c(35, 60)

# Panel label carries the rate, so the number travels with the picture even when the panel is
# lifted out of the figure for a slide.
rate <- setNames(eur$rate, eur$species)
nlib <- setNames(eur$n_libraries, eur$species)
ndet <- setNames(eur$n_detected, eur$species)
# Three short lines beat two long ones: at this panel width the species name plus the common
# name on one line clipped, and a clipped title reads as a rendering fault.
lab  <- sprintf("%s\n%s\n%s of %s libraries  (%s)",
                DRAW, COMMON[DRAW], comma(ndet[DRAW]), comma(nlib[DRAW]),
                percent(rate[DRAW], accuracy = 0.1))
names(lab) <- DRAW

d <- d[d$species %in% DRAW, ]
d$panel <- factor(lab[d$species], levels = lab)

OCEAN <- "white"
world <- map_data("world")
world <- world[world$region != "Antarctica", ]

eu_sites <- sites[sites$lon >= XLIM[1] & sites$lon <= XLIM[2] &
                  sites$lat >= YLIM[1] & sites$lat <= YLIM[2], c("lon", "lat")]

base <- theme_void(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        strip.text       = element_text(face = "bold", size = 11, hjust = 0,
                                        margin = margin(b = 5, t = 7), lineheight = 1.3),
        panel.spacing    = unit(12, "pt"),
        legend.position  = "bottom",
        legend.title     = element_text(size = 11),
        plot.title       = element_text(face = "bold", size = 15, margin = margin(b = 2)),
        plot.subtitle    = element_text(size = 11, colour = "grey30", margin = margin(b = 8)),
        plot.caption     = element_text(size = 9.5, colour = "grey35", hjust = 0,
                                        margin = margin(t = 8)),
        plot.margin      = margin(12, 14, 10, 14))

p <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey92", colour = OCEAN, linewidth = 0.3) +
  # Filled, not open: at this scale an individual site should be followable across panels.
  geom_point(data = eu_sites, aes(lon, lat), inherit.aes = FALSE,
             shape = 21, size = 1.5, stroke = 0.3, colour = "grey45", fill = "grey80") +
  geom_point(data = d, aes(lon, lat, size = n_detected, fill = species),
             shape = 21, stroke = 0.7, colour = "black", alpha = 0.88) +
  scale_fill_manual(values = FILL, guide = "none") +
  scale_size_area(max_size = 7.5, breaks = c(1, 3, 6, 12),
                  name = "libraries with the species detected") +
  coord_fixed(1.45, xlim = XLIM, ylim = YLIM, expand = FALSE) +
  facet_wrap(~ panel, ncol = 3) +
  guides(size = guide_legend(title.position = "top", nrow = 1,
                             override.aes = list(fill = "grey30", colour = "black"))) +
  labs(title = "Where the unreported pathogens are: Europe",
       subtitle = sprintf(
         "%s of the %s mapped wheat libraries, at %s of %s localities.\nGrey points are every sampled locality; coloured points are detections in studies that did not report the species.",
         comma(eur$n_libraries[1]), comma(sum(sites$n_libraries)),
         comma(eur$n_localities[1]), comma(nrow(sites))),
       caption = "The rate in each panel heading is the measured quantity: point density follows sampling effort, so the map shows where a species occurs, not how often.\nAll libraries were collected from visibly diseased plants, so these are not prevalence surfaces.") +
  base

ragg::agg_png(file.path(HERE, "pathogen_map_europe.png"), width = 3900, height = 2750,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "pathogen_map_europe.pdf"), width = 13, height = 9.17)
print(p); invisible(dev.off())
cat("wrote pathogen_map_europe.{png,pdf}\n")
