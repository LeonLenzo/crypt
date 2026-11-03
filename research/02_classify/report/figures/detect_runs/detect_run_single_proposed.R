#!/usr/bin/env Rscript
# detect_run_single_proposed.R — BATCH driver over panel.R.
#
#   Rscript .../detect_run_single_proposed.R [outdir] [none|name|number] [input.tsv]
#
# The drawing moved to panel.R on 2026-10-29 so the Shiny app and this script cannot drift.
# What is left here is the batch loop: read the table, write one PNG and one audit TSV per run
# into runs/<BioProject>/, skipping anything already done.
#
# PREFER THE APP for looking at runs. This exists for the rare case where a whole set has to
# exist as files -- a supplement, or a hand-off. A full cohort render is ~2,862 images, 472 MB
# and 45 minutes, and ANY change to the floors or to taxon_status.tsv invalidates all of it.

args   <- commandArgs(trailingOnly = TRUE)
HERE   <- "02_classify/report/figures/detect_runs"
OUTDIR <- if (length(args) >= 1) args[1] else file.path(HERE, "runs")
LABELS <- if (length(args) >= 2) args[2] else "none"
SRCTSV <- if (length(args) >= 3) args[3] else "detect_runs_proposed.tsv"
stopifnot(LABELS %in% c("name", "none", "number"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(HERE, "panel.R"))

d <- read.delim(file.path(HERE, SRCTSV), stringsAsFactors = FALSE)
d <- d[is.finite(d$redundancy) & is.finite(d$coverage) & d$coverage > 0, ]
BYRUN <- split(seq_len(nrow(d)), d$run)

write_one <- function(run) {
  s <- d[BYRUN[[run]], ]
  p <- make_panel(s, LABELS = LABELS)
  if (is.null(p)) return(invisible(NULL))
  pdir <- file.path(OUTDIR, s$BioProject[1])
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(pdir, paste0(run, ".png"))
  h <- 480 + 52 * attr(p, "n_rows")
  ragg::agg_png(f, width = 3140, height = h, res = 190, background = "white")
  print(p); invisible(dev.off())

  # The audit table: what the panel deliberately drops. Every candidate that cleared the
  # ABSOLUTE floor, not just the calls, because those are the rows where the decision was
  # actually a decision.
  aud <- s[s$call == "yes" | (is.finite(s$distinct) & s$distinct >= 10000), ]
  aud <- aud[order(aud$call != "yes", -aud$coverage), ]
  write.table(
    data.frame(run = aud$run, BioProject = aud$BioProject, host = aud$host,
               species = aud$species, tier = aud$tier, trophic = aud$trophic,
               evidence = aud$evidence, declared = aud$declared, call = aud$call,
               reason = aud$reason, coverage = aud$coverage, saturation = aud$saturation,
               redundancy = aud$redundancy, reads = aud$reads, distinct = aud$distinct),
    file.path(pdir, paste0(run, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("  %-14s %4d candidates, %2d called, %2d undeclared -> %s/%s\n",
              run, nrow(s), attr(p, "n_called"), attr(p, "n_undeclared"),
              s$BioProject[1], basename(f)))
}

runs <- unique(d$run[order(-d$panel_order, d$run)])
if (!exists("FORCE")) FORCE <- nchar(Sys.getenv("RENDER_FORCE")) > 0
proj <- d$BioProject[match(runs, d$run)]
done <- file.exists(file.path(OUTDIR, proj, paste0(runs, ".png"))) &
        file.exists(file.path(OUTDIR, proj, paste0(runs, ".tsv")))
if (!FORCE && any(done)) {
  cat("resuming:", sum(done), "already rendered,", sum(!done), "to go\n")
  runs <- runs[!done]
}
cat("writing", length(runs), "run image(s) to", OUTDIR, "\n")
invisible(lapply(runs, write_one))
write_legend()
