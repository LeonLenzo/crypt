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
xr <- range(d$redundancy); yr <- c(min(d$coverage), 1)
xs <- 10^seq(log10(xr[1]) - 0.3, log10(xr[2]) + 0.3, length.out = 500)
null  <- data.frame(x = xs, y = 1 - exp(-pmin(700, xs)))
shape <- data.frame(x = xs, y = SAT_FLOOR * (1 - exp(-pmin(700, xs))))

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
  cal$n <- NA_integer_
  cal$n[cal$named] <- seq_len(sum(cal$named))
  lab <- cal[cal$named, ]
  other <- cal[!cal$named, ]

  n_und <- sum(cal$declared == "undeclared")
  n_dec <- sum(cal$declared == "declared")
  sub <- sprintf("%s  ·  %s  ·  %s  ·  %s", m$crop, m$host, m$BioProject, m$selection)
  hdr <- sprintf(paste("%d undeclared pathogen%s  ·  %d declared  ·",
                       "%d other fungi above the floor  ·  %d of %d candidates called"),
                 n_und, if (n_und == 1) "" else "s", n_dec,
                 sum(!cal$named & cal$tier != "host"), nrow(cal), nrow(s))

  key <- data.frame(redundancy = xs[1], coverage = yr[1],
                    tier = factor(names(TIER_COLOUR), levels = names(TIER_COLOUR)))

  g <- ggplot(s, aes(redundancy, coverage)) +
    geom_line(data = null,  aes(x, y), inherit.aes = FALSE,
              colour = "grey45", linewidth = 0.8) +
    geom_line(data = shape, aes(x, y), inherit.aes = FALSE,
              colour = INK, linewidth = 0.6, linetype = "longdash") +
    geom_hline(yintercept = COV_FLOOR, linetype = "dashed",
               linewidth = 0.6, colour = INK) +
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
      ncol = 2,
      override.aes = list(size = 3.6, shape = 16,
                          alpha = unname(TIER_ALPHA[names(TIER_LABEL)])))) +
    scale_x_log10(limits = c(xs[1], xs[length(xs)]),
                  breaks = 10^seq(-6, 4, by = 2), labels = seq(-6, 4, by = 2)) +
    scale_y_log10(limits = yr, breaks = 10^seq(-6, 0, by = 2),
                  labels = seq(-6, 0, by = 2)) +
    labs(x = expression(log[10]~"Redundancy"), y = expression(log[10]~"Coverage")) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.grid       = element_blank(),
          axis.line        = element_line(colour = INK, linewidth = 0.7),
          axis.ticks       = element_line(colour = INK, linewidth = 0.7),
          axis.text        = element_text(colour = INK, size = 9),
          axis.title       = element_text(size = 10),
          legend.position  = "bottom",
          legend.text      = element_text(size = 8.4),
          legend.key.size  = unit(10, "pt"),
          legend.margin    = margin(0, 0, 0, 0),
          plot.margin      = margin(4, 6, 4, 4))

  # ── the ranked table: pathogens and host only ───────────────────────────────────────────
  tb <- lab
  tb$y <- -tb$n
  tb$txt <- sprintf("%-2d  %s", tb$n, tb$species)
  tb$cov <- sprintf("%.1f%%", 100 * tb$coverage)
  tb$sat <- sprintf("%.2f", tb$saturation)
  tb$bold <- ifelse(tb$tier == "host", "plain", "bold")

  tg <- ggplot(tb, aes(0, y)) +
    geom_point(aes(x = 0, colour = tier, alpha = tier), shape = 16, size = 3.6) +
    geom_point(data = tb[tb$declared == "declared", ], aes(x = 0), shape = 21,
               size = 6.0, stroke = 0.8, colour = INK, fill = NA) +
    geom_text(aes(x = 0.45, label = txt), hjust = 0, size = 3.1,
              fontface = tb$bold, colour = INK) +
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
    scale_y_continuous(limits = c(-max(STD_ROWS, nrow(tb)) - 0.6, 0.8)) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(4, 4, 4, 2))

  p <- (g | tg) + patchwork::plot_layout(widths = c(1.72, 1)) +
    patchwork::plot_annotation(
      title = run, subtitle = paste0(sub, "\n", hdr),
      caption = paste("ring: pathogen the study declared   ·",
                      "numbered: reportable pathogens and host; other fungi are coloured only\n",
                      "solid curve: null for reads spread uniformly   ·",
                      "long-dash: shape test (saturation 0.05)   ·",
                      "dash: abundance floor (coverage 1%)"),
      theme = theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title      = element_text(face = "bold", size = 17),
        plot.subtitle   = element_text(size = 10.5, colour = "grey25", lineheight = 1.3),
        plot.caption    = element_text(size = 8.8, colour = "grey35", hjust = 0,
                                       lineheight = 1.3),
        plot.margin     = margin(12, 14, 8, 12)))

  f <- file.path(OUTDIR, paste0(run, ".png"))
  h <- 480 + 52 * max(STD_ROWS, nrow(tb))
  ragg::agg_png(f, width = 2300, height = h, res = 190, background = "white")
  print(p); invisible(dev.off())
  cat(sprintf("  %-14s %4d candidates, %2d called, %2d named, %2d undeclared -> %s\n",
              run, nrow(s), nrow(cal), nrow(lab), n_und, basename(f)))
}

runs <- unique(d$run[order(d$panel_order, d$run)])
cat("writing", length(runs), "run image(s) to", OUTDIR, "\n")
invisible(lapply(runs, plot_run))
