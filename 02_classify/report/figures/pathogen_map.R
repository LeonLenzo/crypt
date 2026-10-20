# pathogen_map.R — where the undeclared pathogens are. Styling only.
#
# Every number comes from prep_pathogen_map.py. Six panels, one per species, same basemap and
# same extent so a reader compares positions rather than re-reading axes.
#
# WHICH SIX. The prep qualifies twelve species. Drawn here are the six that are genuine wheat
# LEAF pathogens and together make the classic wheat foliar disease complex: septoria tritici
# blotch, brown rust, powdery mildew, septoria nodorum blotch, snow mould, tan spot. The other
# six stay in pathogen_map.tsv and are named in the caption. Two of those deserve suspicion
# rather than a panel: Cladosporium cladosporioides and Epicoccum nigrum are phyllosphere
# saprotrophs that FUNGuild also labels plant pathogens, and Stemphylium lycopersici is a
# TOMATO pathogen whose presence on wheat is far more likely congeneric cross-mapping from a
# cereal-associated Stemphylium than a real host jump.
#
# BASEMAP follows 01_literature/report/figures/map_crop_origins.R so a map is one map across
# the thesis: maps::map_data("world"), land grey92, and borders drawn in the OCEAN colour so
# they read as hairline separations rather than as drawn lines competing with the points.
#
# THE GREY LAYER IS THE POINT, not decoration. It is every sampled locality, repeated in all
# six panels. Without it an empty region cannot be told from an unsampled one, and for wheat
# the sampling is overwhelmingly northwest European.
#
#   Rscript 02_classify/report/figures/pathogen_map.R

suppressPackageStartupMessages({
  library(ggplot2); library(maps); library(scales); library(ragg)
})

HERE  <- "02_classify/report/figures"
d     <- read.delim(file.path(HERE, "pathogen_map.tsv"), stringsAsFactors = FALSE)
sites <- read.delim(file.path(HERE, "pathogen_map_sites.tsv"), stringsAsFactors = FALSE)

DRAW <- c("Zymoseptoria tritici", "Puccinia triticina", "Blumeria graminis",
          "Parastagonospora nodorum", "Microdochium nivale", "Pyrenophora tritici-repentis")
COMMON <- c("Zymoseptoria tritici"         = "septoria tritici blotch",
            "Puccinia triticina"           = "brown rust",
            "Blumeria graminis"            = "powdery mildew",
            "Parastagonospora nodorum"     = "septoria nodorum blotch",
            "Microdochium nivale"          = "snow mould",
            "Pyrenophora tritici-repentis" = "tan spot")
FILL <- c("Zymoseptoria tritici"         = "#3A6EA5",
          "Puccinia triticina"           = "#A9512C",
          "Blumeria graminis"            = "#45896A",
          "Parastagonospora nodorum"     = "#D9A441",
          "Microdochium nivale"          = "#8C6D8F",
          "Pyrenophora tritici-repentis" = "#8c0404")

d <- d[d$species %in% DRAW, ]
lab <- sprintf("%s\n%s", DRAW, COMMON[DRAW])
names(lab) <- DRAW
d$panel <- factor(lab[d$species], levels = lab)

OCEAN <- "white"
world <- map_data("world")
world <- world[world$region != "Antarctica", ]

# Extent from the DATA, with padding, rather than a whole-globe view: no wheat locality sits
# east of 122 or south of -39, and a full globe would shrink Europe, where the signal is.
XLIM <- c(-128, 130); YLIM <- c(-45, 62)

base <- theme_void(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        strip.text       = element_text(face = "bold", size = 11.5, hjust = 0,
                                        margin = margin(b = 4, t = 6), lineheight = 1.15),
        panel.spacing    = unit(10, "pt"),
        legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.title     = element_text(size = 11),
        plot.title       = element_text(face = "bold", size = 15, margin = margin(b = 2)),
        plot.subtitle    = element_text(size = 11, colour = "grey30", margin = margin(b = 8)),
        plot.caption     = element_text(size = 9.5, colour = "grey35", hjust = 0,
                                        margin = margin(t = 8)),
        plot.margin      = margin(12, 14, 10, 14))

# The denominator layer has no `panel` column, so it is replicated across facets by construction
# rather than by inherit.aes, which would try to map species onto it.
grey_sites <- sites[, c("lon", "lat")]

p <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey92", colour = OCEAN, linewidth = 0.3) +
  geom_point(data = grey_sites, aes(lon, lat), inherit.aes = FALSE,
             shape = 21, size = 1.15, stroke = 0.35, colour = "grey55", fill = NA) +
  geom_point(data = d, aes(lon, lat, size = n_detected, fill = species),
             shape = 21, stroke = 0.65, colour = "black", alpha = 0.85) +
  scale_fill_manual(values = FILL, guide = "none") +
  scale_size_area(max_size = 6.5, breaks = c(1, 3, 6, 12),
                  name = "libraries with the species detected") +
  coord_fixed(1.32, xlim = XLIM, ylim = YLIM, expand = FALSE) +
  facet_wrap(~ panel, ncol = 3) +
  guides(size = guide_legend(title.position = "top", nrow = 1,
                             override.aes = list(fill = "grey30", colour = "black"))) +
  labs(title = "Where the unreported pathogens are",
       subtitle = sprintf(
         "%s wheat libraries with field coordinates, %s localities, %s countries. A point is drawn only where the species was detected in a study that did not report it.",
         comma(sum(sites$n_libraries)), comma(nrow(sites)),
         comma(length(unique(sites$country)))),
       caption = "Open grey points are every sampled locality, repeated in each panel, so an empty region reads as sampled-and-not-found rather than unsampled.\nPoint density follows sampling effort: the map locates a species, it does not measure how often it occurs. For that, see the per-country rates.\nAll libraries were collected from visibly diseased plants, so these are not prevalence surfaces.") +
  base

ragg::agg_png(file.path(HERE, "pathogen_map.png"), width = 4200, height = 2450,
              res = 300, background = "white")
print(p); invisible(dev.off())
cairo_pdf(file.path(HERE, "pathogen_map.pdf"), width = 14, height = 8.17)
print(p); invisible(dev.off())
cat("wrote pathogen_map.{png,pdf}\n")
