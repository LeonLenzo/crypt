#!/usr/bin/env Rscript
# detect_run_single.R — one image per run, for reading the cohort sample by sample.
#
#   Rscript 02_classify/report/figures/detect_run_single.R [outdir]
#
# Reads whatever runs are in detect_runs.tsv and writes <outdir>/<run>.png, one file each.
# Generate the input for any set of runs first:
#   python 02_classify/report/figures/prep_detect_runs.py --runs SRR15141141,SRR8204137
#
# SAME AXES AND SAME BOUNDARIES as detect_runs.R, deliberately. These are meant to be paged
# through next to that figure, so a point cannot mean one thing here and another there.
#
# WHAT IS ADDED: the declared-versus-detected axis. detect_runs.R fills a called point by
# GUILD, which answers "what kind of organism is this" and not "did the originating study say
# it was here" -- and the second is the question the paper is actually asking. A filled dark
# point is a pathogen the study never reported; a ringed one is the pathogen it set out to
# sequence. Reading a run is then: how many solid dark points are there, and how far above the
# floor do they sit.
#
# No new hue is introduced for that. The pathogen red is detect_runs.R's own, and declared
# versus undeclared is carried by FILL versus RING on the same colour, so the two figures stay
# one visual language. Host and non-pathogen fungi are deliberately quiet: they are context.

suppressPackageStartupMessages({ library(ggplot2); library(scales); library(patchwork) })

args   <- commandArgs(trailingOnly = TRUE)
HERE   <- "02_classify/report/figures"
OUTDIR <- if (length(args) >= 1) args[1] else file.path(HERE, "runs")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(file.path(HERE, "detect_runs.tsv"), stringsAsFactors = FALSE)
d <- d[is.finite(d$redundancy) & is.finite(d$coverage) & d$coverage > 0, ]

COV_FLOOR <- 0.01
SAT_FLOOR <- 0.05
INK   <- "#1b1b1b"

# ── HOW THE TAXON ENTERED THE REFERENCE ──────────────────────────────────────────────────────
# Colours follow db_contents.R so the database figure and the per-run figures use one key, and
# PHI-BASE IS KEPT APART FROM FUNGUILD because they are different kinds of evidence: a PHI-base
# record is documented pathogenicity on a named host, a FUNGuild assignment is a guild label.
# db_contents.R splits the target route two ways (recorded on these crops, recorded on another
# plant); the Methods needs a third, the species-rank guild assignment with no host record, and
# that one is given its own HUE rather than a third red so the two sources never read as one
# gradient.
TIER_COLOUR <- c(
  phibase_crop      = "#8c0404",   # db_contents target_cereal
  phibase_plant     = "#e77e7e",   # db_contents target_other
  funguild          = "#D9A441",   # species-rank guild, no PHI-base host record
  host_associate    = "#22127e",   # db_contents
  plant_associate   = "#add0e4",   # db_contents
  target_unresolved = "#bdb5b5",   # target route, not reportable as a pathogen
  congener          = "#dfd8d8",   # db_contents; never callable, kept for completeness
  host              = "#3e9e3b")   # db_contents
TIER_LABEL <- c(
  phibase_crop      = "PHI-base: pathogen of this crop",
  phibase_plant     = "PHI-base: pathogen of another plant",
  funguild          = "FUNGuild: plant pathogen, no host record",
  host_associate    = "host associate: seen in a phytobiome study",
  plant_associate   = "plant associate: recorded on a plant",
  target_unresolved = "pathogen route, not species-resolved",
  congener          = "congener: shares a genus with a pathogen",
  host              = "host (plant)")
# RELEVANCE, carried by alpha instead of by a black ring on every point. The host and a
# pathogen documented on this crop are the two things a reader is here for, so they stay solid;
# everything else fades in order of how much it bears on the question. The ramp is not
# decoration: the bottom two entries are taxa the reference carries to give the classifier
# competitors rather than to be found, and they should recede.
TIER_ALPHA <- c(
  phibase_crop      = 1.00,
  phibase_plant     = 0.85,
  funguild          = 0.70,
  host_associate    = 0.55,
  plant_associate   = 0.42,
  target_unresolved = 0.32,
  congener          = 0.25,
  host              = 1.00)

# Short forms for the table, where the full legend wording would not fit and is redundant with
# the swatch beside it.
TIER_SHORT <- c(
  phibase_crop      = "PHI-base, this crop",
  phibase_plant     = "PHI-base, other plant",
  funguild          = "FUNGuild",
  host_associate    = "host associate",
  plant_associate   = "plant associate",
  target_unresolved = "pathogen route, unresolved",
  congener          = "congener",
  host              = "host")

# Only these are labelled and tabulated. The Methods' definition of a reportable pathogen is
# route 2, or route 3 resolved at or below species rank, which is exactly the first three.
NAMED <- c("phibase_crop", "phibase_plant", "funguild", "host")
# Fixed ROW PITCH for the side table, so a run with 3 calls and a run with 28 are the same
# picture at the same scale rather than two differently-stretched ones. Sized at 30 because the
# cohort's calls per run are median 2, p95 22, p99 31: 99% of images come out identical and
# only the long tail (max 72 fungi in one run) grows the canvas.
STD_ROWS <- 30

# One x range for EVERY run, so paging through is a comparison and not a sequence of
# separately-scaled pictures. Same reason the y range is fixed.
# BOTH AXES ARE FIXED CONSTANTS, Leon's call: x from 10^-6 to 10^4 and y from 10^-6 to 1, NOT
# derived from the data. Deriving them meant the frame moved with whatever happened to be in
# the current detect_runs.tsv, so two runs rendered in different batches were not comparable,
# which is the one thing this figure has to be. A constant frame also puts the null curve, the
# shape test and the abundance floor in the same place on every image, so a reader learns the
# geometry once.
#   Clipping y at -3 was tried and reverted: it removed ~91% of the rejected candidates, and
# that cloud is the evidence that the floors reject far more than they admit.
xr <- c(1e-6, 1e4); yr <- c(1e-6, 1)
xs <- 10^seq(log10(xr[1]), log10(xr[2]), length.out = 500)
# The three boundaries as ONE data frame, so they can be a legend instead of a sentence of
# grey text under the plot. They are distinguished by LINETYPE and share a colour, because
# ggplot allows one colour scale per plot and that one is spent on the tiers.
LINE_COL <- "grey30"
LINE_LAB <- c(null  = "null: reads spread uniformly",
              shape = "shape test (saturation 0.05)",
              floor = "abundance floor (coverage 1%)")
bounds <- rbind(
  data.frame(x = xs, y = 1 - exp(-pmin(700, xs)),                   which = "null"),
  data.frame(x = xs, y = SAT_FLOOR * (1 - exp(-pmin(700, xs))),     which = "shape"),
  data.frame(x = xs, y = COV_FLOOR,                                 which = "floor"))
bounds$which <- factor(bounds$which, levels = names(LINE_LAB))

plot_run <- function(run) {
  s <- d[d$run == run, ]
  if (!nrow(s)) return(invisible(NULL))
  m <- s[1, ]
  rej <- s[s$call != "yes", ]
  cal <- s[s$call == "yes", ]
  cal$tier <- factor(cal$tier, levels = names(TIER_COLOUR))
  cal <- cal[order(-cal$coverage), ]

  # Only reportable pathogens and the host are numbered. The other fungi above the floor are
  # coloured by how they entered the reference and left unlabelled: naming twenty Alternaria
  # and Cladosporium that the Methods does not count as pathogens buries the ones it does.
  cal$named <- cal$tier %in% NAMED
  # NUMBERED IN TIER ORDER, not by coverage alone: the table is read as groups, so the colours
  # cluster and the eye can take "what did PHI-base call on this crop" in one block instead of
  # tracking one hue down an interleaved list. Within a tier it is still coverage, descending.
  cal <- cal[order(cal$tier, -cal$coverage), ]
  cal$n <- NA_integer_
  cal$n[cal$named] <- seq_len(sum(cal$named))
  lab <- cal[cal$named, ]
  other <- cal[!cal$named, ]

  n_und <- sum(cal$declared == "undeclared")
  n_dec <- sum(cal$declared == "declared")
  # Title block, Leon's format: run and project, then the sample's own provenance. The
  # declared pathogen is listed only where the study named one, because "none" is 82% of the
  # cohort and a comma with nothing after it reads as missing data rather than as a finding.
  ttl <- sprintf("%s (%s)", run, m$BioProject)
  bits <- c(m$host, m$declared_pathogens, m$location, m$year, m$tissue)
  sub <- paste(bits[!is.na(bits) & nzchar(bits)], collapse = ", ")
  hdr <- sprintf(paste("%d undeclared pathogen%s  ·  %d declared  ·",
                       "%d other fungi above the floor  ·  %d of %d candidates called"),
                 n_und, if (n_und == 1) "" else "s", n_dec,
                 sum(!cal$named & cal$tier != "host"), nrow(cal), nrow(s))

  key <- data.frame(redundancy = xs[1], coverage = yr[1],
                    tier = factor(names(TIER_COLOUR), levels = names(TIER_COLOUR)))

  g <- ggplot(s, aes(redundancy, coverage)) +
    geom_line(data = bounds, aes(x, y, linetype = which), inherit.aes = FALSE,
              colour = LINE_COL, linewidth = 0.65) +
    scale_linetype_manual(values = c(null = "solid", shape = "longdash", floor = "dashed"),
                          labels = LINE_LAB, breaks = names(LINE_LAB), name = NULL) +
    geom_point(data = rej, shape = 21, size = 1.5, stroke = 0.3,
               colour = "grey60", fill = NA, alpha = 0.5) +
    # EVERY tier needs a row somewhere or its legend key renders empty. drop = FALSE keeps the
    # break in the legend but ggplot draws no glyph for a level no layer supplies, which is how
    # "PHI-base: pathogen of another plant" and "congener" came out as labels with no swatch.
    # An invisible point per tier fixes it, and the key stays complete across every run so the
    # legend does not change shape as you page through.
    geom_point(data = key, aes(colour = tier), shape = 16, size = 0,
               alpha = 0, show.legend = TRUE) +
    # SHAPE 16, not 21. A borderless point has to map COLOUR, because shape 21's colour is its
    # border and setting that to NA with stroke 0 makes ggplot drop the glyph entirely -- the
    # first attempt rendered a plot with numbers and no points at all.
    geom_point(data = other, aes(colour = tier, alpha = tier), shape = 16, size = 3.0) +
    geom_point(data = lab, aes(colour = tier, alpha = tier), shape = 16, size = 4.4) +
    # DECLARED is a halo, not a fill, so it composes with any tier colour instead of
    # competing with it. A pathogen the study named carries a ring; everything else does not.
    geom_point(data = cal[cal$declared == "declared", ], shape = 21, size = 7.0,
               stroke = 0.9, colour = INK, fill = NA) +
    ggrepel::geom_text_repel(data = lab, aes(label = n), size = 3.0, seed = 1,
                             inherit.aes = TRUE, colour = INK,
                             fontface = "bold", max.overlaps = Inf,
                             min.segment.length = 0.1, segment.size = 0.2,
                             segment.colour = "grey70", box.padding = 0.3,
                             point.padding = 0.3) +
    scale_colour_manual(values = TIER_COLOUR, labels = TIER_LABEL,
                        breaks = names(TIER_LABEL), name = NULL, drop = FALSE,
                        aesthetics = "colour") +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none", drop = FALSE) +
    # override.aes must restore alpha: the invisible key layer above is drawn at alpha = 0 and
    # the legend glyph INHERITS that, so any tier with no real data in a run came out as a
    # label with no swatch. Only visible on runs missing a tier, which is most of them.
    # The legend keys carry the SAME alpha as the points, so the ramp is part of the key
    # rather than being flattened to solid and contradicting the plot. override.aes takes a
    # vector, in the order of `breaks`.
    guides(colour = guide_legend(
             ncol = 2, order = 1,
             override.aes = list(size = 3.6, shape = 16,
                                 alpha = unname(TIER_ALPHA[names(TIER_LABEL)]))),
           linetype = guide_legend(ncol = 1, order = 2,
                                   override.aes = list(colour = LINE_COL))) +
    scale_x_log10(limits = c(xs[1], xs[length(xs)]),
                  breaks = 10^seq(-6, 4, by = 2), labels = seq(-6, 4, by = 2)) +
    scale_y_log10(limits = yr, breaks = 10^seq(-6, 0, by = 1),
                  labels = seq(-6, 0, by = 1)) +
    labs(x = expression(log[10]~"Redundancy"), y = expression(log[10]~"Coverage")) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.grid       = element_blank(),
          axis.line        = element_line(colour = INK, linewidth = 0.7),
          axis.ticks       = element_line(colour = INK, linewidth = 0.7),
          axis.text        = element_text(colour = INK, size = 13),
          axis.title       = element_text(colour = INK, size = 15),
          # NO LEGEND ON THE PLOT. Eleven keys in two guide boxes could not be laid out inside
          # the panel width and ggplot clipped the left column. It is also the same legend on
          # every one of ~2,400 images, so it is written ONCE to _legend.png beside them and
          # drawn by hand below, which removes the layout constraint entirely.
          legend.position  = "none",
          legend.margin    = margin(0, 0, 0, 0),
          plot.margin      = margin(4, 6, 4, 4))

  # ── the ranked table, grouped by route ──────────────────────────────────────────────────
  # One block per tier, in the legend's order, with a gap between blocks so the colours read as
  # groups. Tiers that are not individually named still get a row, as a COUNT: a reader must be
  # able to see that eighteen plant associates cleared the floor without being made to read
  # eighteen names the Methods does not count as pathogens.
  rows_tb <- list(); y <- 0
  # `grp`, NOT `sub`: `sub` is the SUBTITLE string, built above and consumed by
  # plot_annotation below, and shadowing it here replaced the whole provenance line with a
  # stray "8". R gives no warning for that, and the plot still renders.
  for (tr in names(TIER_COLOUR)) {
    grp <- cal[cal$tier == tr, ]
    if (!nrow(grp)) next
    if (tr %in% NAMED) {
      for (i in seq_len(nrow(grp))) {
        y <- y - 1
        rows_tb[[length(rows_tb) + 1]] <- data.frame(
          tier = tr, y = y, txt = sprintf("%-2d  %s", grp$n[i], grp$species[i]),
          cov = sprintf("%.1f%%", 100 * grp$coverage[i]),
          sat = sprintf("%.2f", grp$saturation[i]),
          declared = grp$declared[i],
          bold = if (tr == "host") "plain" else "bold", count = FALSE)
      }
    } else {
      y <- y - 1
      rows_tb[[length(rows_tb) + 1]] <- data.frame(
        tier = tr, y = y,
        txt = sprintf("%d × %s", nrow(grp), TIER_SHORT[[tr]]),
        cov = "", sat = "", declared = "", bold = "plain", count = TRUE)
    }
    y <- y - 0.55                      # the break between route blocks
  }
  tb <- do.call(rbind, rows_tb)
  tb$tier <- factor(tb$tier, levels = names(TIER_COLOUR))

  tg <- ggplot(tb, aes(0, y)) +
    geom_point(aes(x = 0, colour = tier, alpha = tier), shape = 16, size = 3.6) +
    geom_point(data = tb[tb$declared == "declared", ], aes(x = 0), shape = 21,
               size = 6.0, stroke = 0.8, colour = INK, fill = NA) +
    geom_text(aes(x = 0.45, label = txt), hjust = 0, size = 3.1,
              fontface = tb$bold,
              colour = ifelse(tb$count, "grey45", INK)) +
    geom_text(aes(x = 8.1, label = cov), hjust = 1, size = 3.0, family = "mono",
              colour = INK) +
    geom_text(aes(x = 9.4, label = sat), hjust = 1, size = 3.0, family = "mono",
              colour = "grey40") +
    annotate("text", x = 8.1, y = 0.1, label = "cov", hjust = 1, size = 2.9,
             colour = "grey40", fontface = "bold") +
    annotate("text", x = 9.4, y = 0.1, label = "sat", hjust = 1, size = 2.9,
             colour = "grey40", fontface = "bold") +
    scale_colour_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none", drop = FALSE) +
    scale_x_continuous(limits = c(-0.5, 9.6)) +
    scale_y_continuous(limits = c(min(-STD_ROWS, min(tb$y)) - 0.8, 0.8)) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(4, 4, 4, 2))

  p <- (g | tg) + patchwork::plot_layout(widths = c(1.72, 1)) +
    patchwork::plot_annotation(
      title = ttl, subtitle = sub,
      theme = theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title      = element_text(face = "bold", size = 17),
        plot.subtitle   = element_text(size = 10.5, colour = "grey25", lineheight = 1.3),
        plot.margin     = margin(12, 14, 8, 12)))

  f <- file.path(OUTDIR, paste0(run, ".png"))
  h <- 480 + 52 * max(STD_ROWS, abs(min(tb$y)))
  ragg::agg_png(f, width = 2300, height = h, res = 190, background = "white")
  print(p); invisible(dev.off())
  cat(sprintf("  %-14s %4d candidates, %2d called, %2d named, %2d undeclared -> %s\n",
              run, nrow(s), nrow(cal), nrow(lab), n_und, basename(f)))
}

# ── the legend, drawn once as its own image ──────────────────────────────────────────────────
# Hand-built rather than extracted from a plot: cowplot's extraction helpers have moved twice
# across versions, and a legend that silently comes back NULL is worse than one drawn
# explicitly. Everything here mirrors the geoms above, so the key cannot drift from the plot.
write_legend <- function() {
  # ONE COLUMN of tiers, lines beside them. Three columns were tried at two different pitches
  # and both overran, because text size is fixed in POINTS while the column offset is in data
  # units: widening the x range to make room shrinks the units and makes the text occupy MORE
  # of them, so the overlap gets worse rather than better. Going vertical removes the
  # constraint instead of negotiating with it.
  tl <- data.frame(tier = factor(names(TIER_LABEL), levels = names(TIER_COLOUR)),
                   lab = unname(TIER_LABEL),
                   x = 0, y = -(seq_along(TIER_LABEL) - 1))
  ll <- data.frame(lab = factor(unname(LINE_LAB), levels = unname(LINE_LAB)),
                   x = 8.0, y = -(seq_along(LINE_LAB) - 1))
  dl <- data.frame(lab = "pathogen the study declared", x = 8.0, y = -4.2)

  g <- ggplot() +
    geom_point(data = tl, aes(x, y, colour = tier, alpha = tier), shape = 16, size = 4.6) +
    geom_text(data = tl, aes(x + 0.22, y, label = lab), hjust = 0, size = 4.0, colour = INK) +
    geom_segment(data = ll, aes(x = x, xend = x + 0.75, y = y, yend = y, linetype = lab),
                 colour = LINE_COL, linewidth = 0.7) +
    geom_text(data = ll, aes(x + 0.95, y, label = lab), hjust = 0, size = 4.0, colour = INK) +
    # The declared halo. It was missing from the key entirely, which left the one ring on the
    # plot unexplained -- and it is the mark that answers the paper's own question, so it is
    # the last thing that should have gone unlabelled.
    geom_point(data = dl, aes(x + 0.37, y), shape = 16, size = 4.6,
               colour = TIER_COLOUR[["phibase_crop"]]) +
    geom_point(data = dl, aes(x + 0.37, y), shape = 21, size = 8.4, stroke = 0.9,
               colour = INK, fill = NA) +
    geom_text(data = dl, aes(x + 0.95, y, label = lab), hjust = 0, size = 4.0, colour = INK) +
    scale_colour_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = c("solid", "longdash", "dashed"), guide = "none") +
    scale_x_continuous(limits = c(-0.3, 15.0), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-7.6, 0.6), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(14, 14, 14, 14))

  f <- file.path(OUTDIR, "_legend.png")
  ragg::agg_png(f, width = 1800, height = 640, res = 190, background = "white")
  print(g); invisible(dev.off())
  cat("  legend -> _legend.png\n")
}

runs <- unique(d$run[order(d$panel_order, d$run)])
cat("writing", length(runs), "run image(s) to", OUTDIR, "\n")
invisible(lapply(runs, plot_run))
write_legend()
