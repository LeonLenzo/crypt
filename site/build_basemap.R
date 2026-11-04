#!/usr/bin/env Rscript
# build_basemap.R -- export the coastline the paper figures already use, for the web map.
#
#   Rscript site/build_basemap.R
#
# The R figures draw `maps::map_data("world")` (map_crop_origins.R:154, geography.R). The web
# map must draw THE SAME outline, or the site and the paper disagree about where a coast is for
# no reason a reader could explain. So this exports that object rather than fetching a basemap:
# one coastline, two renderers, same discipline as panel.R sharing db_contents.R's colours.
#
# Antarctica is dropped, as it is in every figure -- no cereal is grown there and it costs a
# fifth of the vertical extent.
#
# Output is rings of [lon, lat] rounded to 2 dp. At world scale 0.01 deg is ~1 km, which is
# below one screen pixel until zoom 7, and the rounding takes the file from 4.3 MB to ~1 MB.
# Consecutive duplicate points after rounding are dropped, which is most of the saving.

suppressPackageStartupMessages(library(ggplot2))

w <- map_data("world")
w <- w[w$region != "Antarctica", ]

rings <- lapply(split(w, w$group), function(g) {
  x <- round(g$long, 2); y <- round(g$lat, 2)
  keep <- c(TRUE, x[-1] != x[-length(x)] | y[-1] != y[-length(y)])
  if (sum(keep) < 4) return(NULL)          # a ring needs 4 points to be a shape
  paste0("[", paste0("[", x[keep], ",", y[keep], "]", collapse = ","), "]")
})
rings <- Filter(Negate(is.null), rings)

out <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))),
                 "data", "world.json")
if (!length(out) || is.na(out)) out <- "site/data/world.json"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
writeLines(paste0("[", paste(unlist(rings), collapse = ","), "]"), out)

cat(sprintf("world.json  %d rings  %.2f MB\n", length(rings), file.size(out) / 1e6))
