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
PATH  <- "#8c0404"   # as detect_runs.R
QUIET <- "#9FB0BC"   # host and other fungi: present, not the point
INK   <- "#1b1b1b"
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
  cal <- cal[order(-cal$coverage), ]
  cal$n <- seq_len(nrow(cal))
  cal$kind <- factor(cal$declared,
                     levels = c("undeclared", "declared", "other fungus", "host"))

  # NUMBERS, NOT NAMES, on the plot. A run with 28 calls cannot carry 28 species labels: the
  # first attempt collided Cercospora over Pseudopithomyces and stacked three Bipolaris on the
  # 1% line. A one- or two-character label repels cleanly at any density, and the names move to
  # a ranked table beside the plot where they are easier to read anyway. Numbering is by
  # coverage, descending, so #1 is always the most abundant thing called.
  n_und <- sum(cal$declared == "undeclared")
  n_dec <- sum(cal$declared == "declared")
  sub <- sprintf("%s  ·  %s  ·  %s  ·  %s", m$crop, m$host, m$BioProject, m$selection)
  hdr <- sprintf("%d undeclared pathogen%s  ·  %d declared  ·  %d of %d candidates called",
                 n_und, if (n_und == 1) "" else "s", n_dec, nrow(cal), nrow(s))

  pt <- function(k, ...) cal[cal$kind == k, ]

  g <- ggplot(s, aes(redundancy, coverage)) +
    geom_line(data = null,  aes(x, y), inherit.aes = FALSE,
              colour = "grey45", linewidth = 0.8) +
    geom_line(data = shape, aes(x, y), inherit.aes = FALSE,
              colour = INK, linewidth = 0.6, linetype = "longdash") +
    geom_hline(yintercept = COV_FLOOR, linetype = "dashed",
               linewidth = 0.6, colour = INK) +
    geom_point(data = rej, shape = 21, size = 1.5, stroke = 0.3,
               colour = "grey60", fill = NA, alpha = 0.5) +
    geom_point(data = pt("declared"), shape = 21, size = 3.8,
               stroke = 1.4, colour = PATH, fill = "white") +
    geom_point(data = pt("host"), shape = 21, size = 3.0,
               stroke = 0.5, colour = INK, fill = QUIET) +
    geom_point(data = pt("other fungus"), shape = 21, size = 3.0,
               stroke = 0.5, colour = INK, fill = "white") +
    geom_point(data = pt("undeclared"), shape = 21, size = 3.8,
               stroke = 0.6, colour = INK, fill = PATH) +
    ggrepel::geom_text_repel(data = cal, aes(label = n), size = 3.0, seed = 1,
                             colour = INK, fontface = "bold", max.overlaps = Inf,
                             min.segment.length = 0.1, segment.size = 0.2,
                             segment.colour = "grey70", box.padding = 0.3,
                             point.padding = 0.25) +
    scale_x_log10(limits = c(xs[1], xs[length(xs)]),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(limits = yr, breaks = 10^(-6:0),
                  labels = c("0.0001%", "0.001%", "0.01%", "0.1%", "1%", "10%", "100%")) +
    annotation_logticks(sides = "bl", size = 0.3, short = unit(2, "pt"),
                        mid = unit(3, "pt"), long = unit(4, "pt")) +
    labs(x = "Redundancy available  (reads × k-mers per spot / taxon's db_v4 minimizers)",
         y = "Coverage observed  (distinct / taxon's db_v4 minimizers)") +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.grid       = element_blank(),
          axis.line        = element_line(colour = INK, linewidth = 0.7),
          axis.ticks       = element_line(colour = INK, linewidth = 0.7),
          axis.text        = element_text(colour = INK, size = 9),
          axis.title       = element_text(size = 10),
          plot.margin      = margin(4, 6, 4, 4))

  # ── the ranked table ────────────────────────────────────────────────────────────────────
  # Same glyphs as the plot, so the two halves are read with one key.
  tb <- cal
  tb$y <- -tb$n
  tb$lab <- sprintf("%-2d  %s", tb$n, tb$species)

  tb$cov <- sprintf("%.1f%%", 100 * tb$coverage)
  tb$sat <- sprintf("%.2f", tb$saturation)

  tp <- function(k) tb[tb$kind == k, ]
  tg <- ggplot(tb, aes(0, y)) +
    geom_point(data = tp("declared"), aes(x = 0), shape = 21, size = 3.4,
               stroke = 1.3, colour = PATH, fill = "white") +
    geom_point(data = tp("host"), aes(x = 0), shape = 21, size = 2.7,
               stroke = 0.5, colour = INK, fill = QUIET) +
    geom_point(data = tp("other fungus"), aes(x = 0), shape = 21, size = 2.7,
               stroke = 0.5, colour = INK, fill = "white") +
    geom_point(data = tp("undeclared"), aes(x = 0), shape = 21, size = 3.4,
               stroke = 0.6, colour = INK, fill = PATH) +
    geom_text(aes(x = 0.35, label = lab), hjust = 0, size = 3.1,
              fontface = ifelse(tb$kind %in% c("undeclared", "declared"), "bold", "plain"),
              colour = ifelse(tb$kind %in% c("undeclared", "declared"), PATH, INK)) +
    geom_text(aes(x = 8.1, label = cov), hjust = 1, size = 3.0, family = "mono",
              colour = INK) +
    geom_text(aes(x = 9.4, label = sat), hjust = 1, size = 3.0, family = "mono",
              colour = "grey40") +
    annotate("text", x = 8.1, y = 0.1, label = "cov", hjust = 1, size = 2.9,
             colour = "grey40", fontface = "bold") +
    annotate("text", x = 9.4, y = 0.1, label = "sat", hjust = 1, size = 2.9,
             colour = "grey40", fontface = "bold") +
    scale_x_continuous(limits = c(-0.4, 9.6)) +
    scale_y_continuous(limits = c(-max(STD_ROWS, nrow(tb)) - 0.6, 0.8)) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(4, 4, 4, 2))

  p <- (g | tg) + patchwork::plot_layout(widths = c(1.72, 1)) +
    patchwork::plot_annotation(
      title = run, subtitle = paste0(sub, "\n", hdr),
      caption = paste("filled red: pathogen the study did not report   ·",
                      "ringed red: pathogen it declared   ·",
                      "grey fill: host   ·  white fill: other fungus\n",
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
  h <- 420 + 52 * max(STD_ROWS, nrow(cal))
  ragg::agg_png(f, width = 2300, height = h, res = 190, background = "white")
  print(p); invisible(dev.off())
  cat(sprintf("  %-14s %4d candidates, %2d called, %2d undeclared -> %s\n",
              run, nrow(s), nrow(cal), n_und, basename(f)))
}

runs <- unique(d$run[order(d$panel_order, d$run)])
cat("writing", length(runs), "run image(s) to", OUTDIR, "\n")
invisible(lapply(runs, plot_run))
