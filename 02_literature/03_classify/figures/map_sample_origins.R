#!/usr/bin/env Rscript
# map_sample_origins.R — where the analysed cohort was collected, and when.
#
# Proportional-symbol world map: one circle per country, area = number of BioSamples,
# fill = the period its median collection year falls in. A companion panel gives the
# per-year distribution, because one fill per country cannot show that the UK spans
# 1978 to 2022.
#
# Area, not radius, encodes count (scale_size_area), so a circle twice the area means
# twice the samples. Radius-scaling would overstate the large countries ~2x.
#
# Period is an ORDERED band, so fill is a single-hue sequential ramp light to dark, not a
# categorical palette and not a rainbow. Countries whose samples carry no collection date
# are grey rather than dropped, so the map does not silently understate them.
#
# Circles sit at the centroid of each country's LARGEST polygon, not the midpoint of its
# bounding box. The bounding box is wrong for any country with distant territory: it puts the
# USA in the Pacific (Alaska to Maine) and France in the Atlantic (metropolitan plus French
# Guiana and Reunion).
#
# Country-level origin only. Within-country structure exists for the rust surveys (PRJEB39201
# has 235 localities, PRJEB31334 75) but only ~7% of the cohort carries lat/lon, so plotting
# true points would misrepresent coverage.
#
# Exact counts live in the ranked bar panel rather than as labels on the map: 20 of the 53
# countries are European and any on-map labelling collides there.
#
#   Rscript 02_literature/03_classify/figures/map_sample_origins.R
# Output: figures/map_sample_origins.{png,pdf}

suppressMessages({
  library(ggplot2); library(maps); library(dplyr); library(scales)
  library(patchwork); library(ragg)
})

here  <- "02_literature/03_classify"
dat   <- read.delim(file.path(here, "data/country_summary.tsv"), stringsAsFactors = FALSE)
prov  <- read.delim(file.path(here, "data/cohort_provenance.tsv"), stringsAsFactors = FALSE)

# our normalised names -> the names the maps package uses for its regions
to_maps <- c("United Kingdom" = "UK", "Czechia" = "Czech Republic")
dat$region <- ifelse(dat$country %in% names(to_maps), to_maps[dat$country], dat$country)

world <- map_data("world")
world <- world[world$region != "Antarctica", ]

# centroid of the largest polygon per country, so outlying territory cannot drag the point
cent <- world %>%
  group_by(region, group) %>%
  summarise(long = mean(long), lat = mean(lat), npts = n(), .groups = "drop") %>%
  group_by(region) %>%
  slice_max(npts, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(region, long, lat)

pts <- dat %>%
  left_join(cent, by = "region") %>%
  filter(!is.na(long))

dropped <- dat %>% anti_join(cent, by = "region")
if (nrow(dropped)) message("unplaced: ", paste(dropped$country, collapse = ", "))

PERIODS <- c("pre-2000", "2000-2009", "2010-2014", "2015-2019", "2020-2025", "unknown")
pts$period <- factor(pts$period, levels = PERIODS)

# single-hue sequential blues, light to dark; grey sits outside the ramp for "no date"
fills <- c("pre-2000"  = "#c6dbef", "2000-2009" = "#9ecae1", "2010-2014" = "#6baed6",
           "2015-2019" = "#3182bd", "2020-2025" = "#08519c", "unknown"   = "#bdbdbd")

base <- theme_minimal(base_size = 16) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

p_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey94", colour = "grey80", linewidth = 0.3) +
  geom_point(data = pts, aes(long, lat, size = n, fill = period),
             shape = 21, colour = "black", stroke = 1, alpha = 0.85) +
  scale_size_area(max_size = 22, breaks = c(10, 50, 150, 300, 460),
                  labels = comma, name = "BioSamples") +
  # drop unused bands: no country has a pre-2000 MEDIAN (76 samples are pre-2005 but they are
  # spread thinly), and an empty legend key reads as a rendering fault
  scale_fill_manual(values = fills, name = "Median collection year", drop = TRUE) +
  coord_fixed(1.35, xlim = c(-165, 180), ylim = c(-56, 78), expand = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 6), order = 2),
         size = guide_legend(order = 1)) +
  labs(title = "Where the analysed cohort was collected",
       subtitle = sprintf("%s of %s BioSamples with a resolved country, across %d countries",
                          comma(sum(pts$n)), comma(nrow(prov)), nrow(pts))) +
  base +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), legend.position = "right",
        legend.box = "vertical")

# ranked bars carry the exact counts the map cannot, using the same fill scale
top <- pts %>% arrange(desc(n)) %>% head(16) %>%
  mutate(country = factor(country, levels = rev(country)))

p_bar <- ggplot(top, aes(n, country, fill = period)) +
  geom_col(colour = "black", linewidth = 0.6, width = 0.75) +
  geom_text(aes(label = comma(n)), hjust = -0.25, size = 4.4, colour = "grey20") +
  scale_fill_manual(values = fills, guide = "none", drop = TRUE) +
  scale_x_continuous(expand = expansion(c(0, 0.16))) +
  labs(x = "BioSamples", y = NULL, subtitle = "Top 16 of 53 countries") +
  base + theme(panel.grid.major.y = element_blank(),
               panel.grid.minor = element_blank())

# the long pre-2000 tail is 28 samples; cutting the axis there would hide them, so it is
# annotated instead and the axis starts where the data actually begins
yr <- prov %>% filter(!is.na(collection_year), collection_year != "") %>%
  mutate(y = as.integer(collection_year)) %>% count(y)
early <- sum(yr$n[yr$y < 2005])

p_yr <- ggplot(yr %>% filter(y >= 2005), aes(y, n)) +
  geom_col(fill = "#3182bd", colour = "black", linewidth = 0.5, width = 0.8) +
  scale_x_continuous(breaks = seq(2005, 2025, 5)) +
  scale_y_continuous(expand = expansion(c(0, 0.12))) +
  labs(x = "Collection year", y = "BioSamples",
       subtitle = sprintf("%s samples dated 2005 onward; %d earlier (to 1925); %s undated",
                          comma(sum(yr$n) - early), early,
                          comma(nrow(prov) - sum(yr$n)))) +
  base + theme(panel.grid.minor = element_blank(),
               panel.grid.major.x = element_blank())

out <- p_map / (p_bar | p_yr) + plot_layout(heights = c(2.2, 1.3))

agg_png(file.path(here, "figures/map_sample_origins.png"), width = 15, height = 12,
        units = "in", res = 300, background = "white")
print(out); invisible(dev.off())
cairo_pdf(file.path(here, "figures/map_sample_origins.pdf"), width = 15, height = 12)
print(out); invisible(dev.off())
cat("wrote figures/map_sample_origins.{png,pdf}\n")
cat(sprintf("  %d countries plotted, %d samples\n", nrow(pts), sum(pts$n)))
