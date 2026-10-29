# undeclared_map.R -- how many pathogens a library carries that its own study did not report,
# one image per crop, one panel per continent. Styling only; numbers come from
# prep_undeclared_map.py.
#
# ONE IMAGE PER CROP. The crops are not equally far through classification, and putting them on
# one map invites a reader to compare a nearly-finished host with a barely-started one. Wheat is
# effectively done; maize, sorghum and rice are not. Leon's call, 2026-10-26.
#
# WHOLE CONTINENTS, NOT BOXES AROUND THE POINTS. Cropping to the data makes sparse sampling look
# like thorough sampling of a small place, and the unsampled interior is exactly what a reader
# needs to see.
#
# WHY COLOUR IS A PER-LIBRARY MEDIAN AND SIZE IS EFFORT. A locality with forty libraries finds
# more distinct species than one with two for reasons of effort alone. Effort already sets point
# SIZE, so colouring by the distinct-species count would encode the same confound twice.
#
# ZEROS ARE A DIFFERENT MARK, NOT THE BOTTOM OF THE RAMP. 40.5% of wheat localities found
# nothing unreported and a further 33.8% sit in the lowest bin, so with a solid-dark low end
# three quarters of the map was dark and two different statements looked identical. Zeros are
# now hollow grey rings and the ramp starts at viridis 0.18, well clear of near-black, so "one
# per library" reads as a colour rather than as absence. Leon read the old version as zeros,
# which is the figure's fault.
#
# NO BLACK OUTLINES ON THE DATA POINTS. At this density a black stroke on every point reads as
# crowding rather than as definition. The outline is reserved for the zero ring, where it is
# the whole signal. Fills are alpha 0.6 so overlapping localities show as density.
#
#   Rscript 02_classify/report/figures/undeclared_map/undeclared_map.R

suppressPackageStartupMessages({
  library(ggplot2); library(maps); library(scales); library(ragg); library(patchwork)
})

HERE <- "02_classify/report/figures/undeclared_map"
d    <- read.delim(file.path(HERE, "undeclared_map.tsv"), stringsAsFactors = FALSE)
reg  <- read.delim(file.path(HERE, "undeclared_map_regions.tsv"), stringsAsFactors = FALSE)

# ── GLOBAL SCALES ───────────────────────────────────────────────────────────────────────────
# Every panel is also written as its own image, so the scales must be set ONCE here and shared
# by all of them. A per-image scale would make two images that look alike mean different
# things, which is the specific failure separate images invite.
#
# COLOUR spans every crop, so a maize panel and a wheat panel are directly comparable. Bins are
# absolute, so wheat simply leaves the top bins unused rather than being rescaled into them.
CAP    <- max(d$median_undeclared)                  # 18, from maize at Pendleton
CBREAK <- c(1, 2, 3, 4, 6, 9)
#
# SIZE is capped rather than spanning the full range. One sorghum locality holds 284 libraries
# against a wheat maximum of 30 and a wheat 90th percentile of 6, so a scale reaching 284 would
# draw almost every other locality as a dot. The cap is 60 and anything above it is drawn at
# the maximum size; exactly one locality in the cohort is affected and the caption says so.
SCAP   <- 60
SBREAK <- c(1, 5, 20, 60)

OCEAN <- "white"
LAND  <- "grey85"
MIN_LOCALITIES <- 2          # below this a map is a list of two dots; the count is printed
CAP_TXT <- paste("Hollow grey rings are localities where nothing unreported was found;",
                 "filled points carry at least one.\nPoint size is sampling effort, not signal.",
                 "Colour and size scales are identical across every crop and continent.")

world <- map_data("world")
world <- world[world$region != "Antarctica", ]

base <- theme_void(base_size = 13) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        plot.title       = element_text(face = "bold", size = 12.5,
                                        margin = margin(b = 2, t = 2)),
        legend.position  = "bottom",
        legend.title     = element_text(size = 11),
        plot.margin      = margin(4, 6, 4, 6))

panel <- function(dd, xlim, ylim, title, ratio = 1.3) {
  zero <- dd[dd$median_undeclared == 0, ]
  hit  <- dd[dd$median_undeclared >  0, ]
  ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = LAND, colour = OCEAN, linewidth = 0.3) +
    geom_point(data = zero, aes(lon, lat), inherit.aes = FALSE,
               shape = 21, size = 1.7, stroke = 0.55, colour = "grey45", fill = NA,
               alpha = 0.9) +
    geom_point(data = hit, aes(lon, lat, colour = median_undeclared, size = n_libraries),
               shape = 16, alpha = 0.6) +
    # BINNED, not continuous. Most localities sit at one or two per library against a cap of
    # eight, so a linear ramp spends most of its range where there is no data and everything
    # reads as the same dark purple. Bins put the contrast where the values actually are.
    scale_colour_viridis_b(limits = c(0, CAP), option = "viridis", begin = 0.18,
                           breaks = CBREAK, oob = scales::squish,
                           name = "unreported pathogens per library (median)") +
    # limits, not just breaks: scale_size_area drops breaks outside a panel's own range, which
    # gives every panel a different scale and leaves patchwork unable to merge the legends.
    scale_size_area(max_size = 7.5, limits = c(0, SCAP), breaks = SBREAK,
                    oob = scales::squish, name = "libraries at the locality") +
    coord_fixed(ratio, xlim = xlim, ylim = ylim, expand = FALSE) +
    guides(colour = guide_coloursteps(title.position = "top", barwidth = 13, barheight = 0.6,
                                      order = 1, show.limits = TRUE),
           size   = guide_legend(title.position = "top", nrow = 1, order = 2,
                                 override.aes = list(colour = "grey40", alpha = 1))) +
    labs(title = title) +
    base
}

# Rows of at most three, each row's widths proportional to panel aspect so continents of very
# different shape still tile evenly.
rows_of <- function(pan, wid, ncol = 3) {
  idx <- split(seq_along(pan), ceiling(seq_along(pan) / ncol))
  lapply(idx, function(i) wrap_plots(pan[i], nrow = 1, widths = wid[i]))
}

for (crop in sort(unique(d$crop))) {
  dd <- d[d$crop == crop, ]
  rr <- reg[reg$crop == crop, ]
  if (nrow(dd) < MIN_LOCALITIES || !nrow(rr)) {
    cat(sprintf("  %-14s %d localities, not drawn\n", crop, nrow(dd))); next
  }
  if (max(dd$median_undeclared) == 0) {
    cat(sprintf("  %-14s nothing undeclared, not drawn\n", crop)); next
  }

  pan <- lapply(seq_len(nrow(rr)), function(i) {
    r <- rr[i, ]
    sub <- dd[dd$lon >= r$xmin & dd$lon <= r$xmax & dd$lat >= r$ymin & dd$lat <= r$ymax, ]
    panel(sub, c(r$xmin, r$xmax), c(r$ymin, r$ymax), r$region)
  })
  wid <- ((rr$xmax - rr$xmin) * 1.3) / (rr$ymax - rr$ymin)

  rws <- rows_of(pan, wid)
  out <- Reduce(`/`, rws) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s: pathogens the study did not report", crop),
      caption = "Hollow grey rings are localities where nothing unreported was found; filled points carry at least one. Point size is sampling effort, not signal.\nColour and size scales are identical across every crop and continent. One locality holds more than 60 libraries and is drawn at the maximum point size.",
      theme = theme(plot.title   = element_text(face = "bold", size = 17,
                                                margin = margin(b = 6)),
                    plot.caption = element_text(size = 9.5, colour = "grey35", hjust = 0,
                                                margin = margin(t = 6)),
                    plot.background = element_rect(fill = "white", colour = NA),
                    legend.position = "bottom", legend.box = "horizontal"))

  # Each continent also as its own standalone image, on the same global scales.
  slug <- tolower(gsub(" ", "_", crop))
  for (i in seq_len(nrow(rr))) {
    r <- rr[i, ]
    sub <- dd[dd$lon >= r$xmin & dd$lon <= r$xmax & dd$lat >= r$ymin & dd$lat <= r$ymax, ]
    one <- panel(sub, c(r$xmin, r$xmax), c(r$ymin, r$ymax), r$region) +
      labs(title = sprintf("%s, %s: pathogens the study did not report", crop, r$region),
           caption = CAP_TXT) +
      theme(plot.title   = element_text(face = "bold", size = 14, margin = margin(b = 6)),
            plot.caption = element_text(size = 9, colour = "grey35", hjust = 0,
                                        margin = margin(t = 6)))
    asp <- ((r$ymax - r$ymin)) / ((r$xmax - r$xmin) * 1.3)
    w <- 2800; h <- round(w * asp) + 950
    rslug <- tolower(gsub(" ", "_", r$region))
    ragg::agg_png(file.path(HERE, sprintf("undeclared_map_%s_%s.png", slug, rslug)),
                  width = w, height = h, res = 300, background = "white")
    print(one); invisible(dev.off())
    cairo_pdf(file.path(HERE, sprintf("undeclared_map_%s_%s.pdf", slug, rslug)),
              width = w / 300, height = h / 300)
    print(one); invisible(dev.off())
  }

  h <- 1200 + 1100 * length(rws)
  ragg::agg_png(file.path(HERE, sprintf("undeclared_map_%s.png", slug)),
                width = 4000, height = h, res = 300, background = "white")
  print(out); invisible(dev.off())
  cairo_pdf(file.path(HERE, sprintf("undeclared_map_%s.pdf", slug)),
            width = 13.3, height = h / 300)
  print(out); invisible(dev.off())
  cat(sprintf("  %-14s %d localities, %d libraries, %d panels -> undeclared_map_%s.png\n",
              crop, nrow(dd), sum(dd$n_libraries), nrow(rr), slug))
}
