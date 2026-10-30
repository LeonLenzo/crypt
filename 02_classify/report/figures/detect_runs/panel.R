# panel.R — the per-run detection panel, as a FUNCTION rather than a script.
#
# Extracted 2026-10-29. The batch renderer drew 2,862 PNGs in two hours and a 23-species
# change to taxon_status.tsv invalidated all of them the same afternoon. Drawing on demand
# from the data is the right shape; drawing thousands of files each time a threshold moves is
# not. So the drawing lives here and has two consumers:
#
#   detect_run_single_proposed.R   batch, writes PNG + audit TSV per run
#   app/app.R                      Shiny, renders one run live with the floors as sliders
#
# Nothing in this file reads a table or writes a file. `make_panel()` takes the candidate rows
# for ONE run and returns a patchwork object; the caller decides what to do with it.
#
# FLOORS ARE ARGUMENTS, not constants, which is what lets the app move them. The defaults are
# 05_detect/detect.py's own and must track it: COV_FLOOR 0.02, SAT_FLOOR 0.05, MIN_DISTINCT
# 10,000. If detect.py changes, change them here.

suppressPackageStartupMessages({ library(ggplot2); library(scales); library(patchwork) })


COV_FLOOR <- 0.02
SAT_FLOOR <- 0.05
INK   <- "#1b1b1b"

# ── WHAT THE ORGANISM IS TO THIS PLANT ───────────────────────────────────────────────────────
# Two questions, answered in order, and the order is the whole design. FIRST: is it a pathogen
# at all? That is species-rank PHI-base, or species-rank FUNGuild, or the curated verdict in
# taxon_status.tsv, which overrides both where it reaches. SECOND, and independently: has it
# ever been recorded on the crop in front of us? That is USDA `recorded_on`, matched at genus.
#
# Crossing them gives four fungal keys plus grey for no evidence either way, plus the host.
# Hue carries the first question, depth within the hue carries the second, so a reader learns
# "red is a pathogen, dark red is a pathogen OF THIS CROP" once and the key stops being eight
# things to memorise.
#
# HOST FIRST in every vector, Leon 2026-10-28. It is the one row on the page that is not a
# finding: it says which plant this library came from, and a reader checks it before reading
# anything else. Order is load-bearing here -- factor levels, the legend breaks and the
# override.aes vectors all index these names -- so the five vectors move together.
#
# Colours are db_contents.R's own, unchanged, so the database figure and the per-run figures
# still share one key.
TIER_COLOUR <- c(
  host           = "#3e9e3b",   # db_contents
  pathogen_crop  = "#8c0404",   # db_contents target_cereal
  pathogen_other = "#e77e7e",   # db_contents target_other
  resident_crop  = "#22127e",   # db_contents host_associate
  resident_other = "#add0e4",   # db_contents plant_associate
  no_record      = "#bdb5b5")   # no plant host record; evidence absent, not negative
TIER_LABEL <- c(
  host           = "host (plant)",
  pathogen_crop  = "host pathogen",
  pathogen_other = "non-host pathogen",
  resident_crop  = "resident, recorded on this crop",
  resident_other = "resident, recorded on other plants",
  no_record      = "resident, no plant host record")

# SECTION HEADINGS for the side table: sentence case, bold, and carrying the swatch. The rows
# under a heading no longer repeat it, so a colour is stated once per block instead of once
# per line.
#   HOST / HOST PATHOGEN / NON-HOST PATHOGEN / RESIDENT, Leon 2026-10-28. Four words, and the
# host axis is now in the NAME of the category rather than in a phrase after it, so the table
# of contents for the figure is the argument the figure makes: what the plant is, what
# attacks it, what attacks something else, and what merely lives there.
TIER_HEAD <- c(
  host           = "Host",
  pathogen_crop  = "Host pathogen",
  pathogen_other = "Non-host pathogen",
  resident_crop  = "Recorded on this crop",
  resident_other = "Recorded on other plants",
  no_record      = "No plant host record")

# FILLED versus HOLLOW, a second channel that costs no hue. Leon, 2026-10-28. Filled means the
# taxon is called a pathogen, or is the host; hollow means it is not. Hue then says which
# plant it belongs to and fill says what it is doing there, so the two questions the tier
# crosses are carried by two marks instead of by six colours.
#   `no_record` is hollow with the residents rather than filled. It is the weaker claim of the
# two, not the stronger one, and a filled grey would put it on the pathogen side of the only
# rule the reader has to learn.
HOLLOW <- c("resident_crop", "resident_other", "no_record")
# ONE HEADING over all three hollow tiers, with the colours moved down onto the rows beneath
# it. Leon, 2026-10-28. Three separate single-line blocks made three sections out of one idea
# and gave the non-findings as much structural weight as the findings.
   # "RESIDENT", Leon's call 2026-10-28, chosen over "Not a pathogen". Worth knowing what it
# costs so nobody re-opens it blind: the third row under it, `no_record`, is an ABSENCE of any
# host record rather than a positive record of residency, so the heading asserts for that row
# slightly more than the evidence holds. The row itself still says exactly what is known, and
# the gain is that all four headings now name a kind of organism rather than a kind of
# database entry.
HOLLOW_HEAD <- "Resident"
# Lowercase forms, because these now sit as rows under a heading rather than as headings.
TIER_ROW <- c(
  resident_crop  = "recorded on this crop",
  resident_other = "recorded on other plants",
  no_record      = "no plant host record")
# ── SHAPE = WHAT KIND OF ORGANISM ────────────────────────────────────────────────────────────
# Two independent questions, two independent channels. Colour answers what the taxon is to THIS
# plant; shape answers what kind of thing it is. Leon, 2026-10-28, and shape is the right
# channel precisely because the two do not interact: a biotroph is a biotroph whether it sits
# on its host or on the wrong one, so the figure can now say both at once.
#
# The split that matters is OBLIGATE BIOTROPHY. Rusts, powdery mildews and smuts have no
# saprotrophic phase, so in an RNA-seq library their transcripts require living infected
# tissue. That is the one class of detection the standing caveat -- DNA presence is not proof
# of active infection, because spore drift and a mis-called host also deposit real DNA -- does
# not apply to. 51.3% of host-pathogen calls, 492 of them undeclared across 440 runs.
#   Rejected on the way here: FungalTraits `Plant_pathogenic_capacity` (3,115 of 3,150 calls
# are leaf/fruit/seed_pathogen, because the cohort is foliar, so it splits nothing) and a
# major/minor severity axis (no machine-readable source; EPPO is a regulatory register and
# would rank quarantine status, not importance).
#
# Filled still means pathogen and hollow still means not, so the earlier rule survives: every
# pathogen shape is solid and the resident is the only open one.
# FOUR pathogen classes, ordered by how much is known about the species, Leon 2026-10-29.
# The middle used to split on FungalTraits' genus lifestyle, which put Epicoccum nigrum beside
# Zymoseptoria tritici with only an A/B/C code between them. It now splits on species-rank
# evidence -- PHI-base or a cited curation -- and the two uncharacterised classes are kept
# apart by whether the GENUS is a pathogen genus at all.
#   Shape 25 is a FILLED triangle-down and takes its interior from `fill`, not `colour`, so
# the geoms map both. Residents use shape 1, an open circle that needs no fill, which is what
# frees fill for 25. Swapping 21 for 1 looks identical and removes the clash.
MARK_SHAPE <- c(biotroph    = 17L,  # triangle up: cannot be dead, so detection = infection
                necrotroph  = 16L,  # circle: necrotroph/hemibiotroph with published pathology
                weak        = 15L,  # square: same lifestyle, nothing known at species rank
                opportunist = 25L,  # triangle down: not even a plant-pathogen genus
                host        = 18L,  # diamond; not a fungus, so not a fungal shape
                resident    =  1L)  # the only open one
MARK_LABEL <- c(biotroph    = "obligate biotroph",
                necrotroph  = "necrotroph / hemibiotroph",
                weak        = "weak pathogen (no species-rank record)",
                opportunist = "opportunist (saprotroph genus)",
                host        = "host plant",
                resident    = "not a pathogen")
TROPHIC_ORDER <- c("biotroph", "necrotroph", "weak", "opportunist")

# NOTHING IS HIDDEN, same as before and for the same reason: the grey bucket is 15% of called
# rows and dropping it from a picture whose job is to show what cleared the floors would be a
# lie by omission. What changed is its NAME. It used to read "pathogen genus label only",
# which asserts a genus-level pathogen label; the honest statement is that the reference
# databases hold no plant host record for the species, which is a fact about the databases and
# not about the biology.
HIDE <- character(0)

# ONE ALPHA FOR EVERY TIER, Leon 2026-10-28. The ramp is gone. It was a third channel trying
# to say what hue and fill now say between them -- which plant, and pathogen or not -- and a
# point's importance was being encoded twice over while its opacity also had to be read as
# overplotting. A flat 0.7 means density is the only thing transparency reports, so a stack of
# points in the cloud reads as a stack rather than as a tier.
#   Kept as a named vector rather than a bare constant because the colour guide's override.aes
# indexes it by tier name, and because this is the obvious knob to turn next.
FLAT_ALPHA <- 0.7
TIER_ALPHA <- setNames(rep(FLAT_ALPHA, length(TIER_COLOUR)), names(TIER_COLOUR))

# NAMED, and this is the substantive change to what gets read by name. The old rule named
# every reportable pathogen, where "reportable" meant PHI-base or species-rank FUNGuild. The
# new rule names every PATHOGEN, both halves of the host axis, because a pathogen of another
# host is the drift explanation and a reader cannot weigh it without the species name. The
# non-pathogen tiers stay counted, not listed, exactly as before -- but they now get a heading
# of their own rather than a bare tally line.
NAMED <- c("host", "pathogen_crop", "pathogen_other")
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
#   Y RUNS TO -5, NOT -6. Leon, 2026-10-28. Clipping at -3 was tried and reverted because it
# removed ~91% of the rejected candidates, and that cloud is the evidence that the floors
# reject far more than they admit -- but -5 is not the same trade. Measured over the six
# example runs it hides 6 of 2,551 rejected points (0.2%) and 0 called ones, so the cloud
# survives intact and the decade it buys goes to the band between the floors and the null,
# which is where every decision on the figure is actually made.
#   Re-measure before moving it again. The number that matters is what fraction of the
# REJECTED cloud disappears, not how empty the bottom of the panel looks.
# TEXT SCALE. Every type size below is a base size times this, so the whole image scales with
# one number. 1.3 was set after the canvas widened from 2300 to 2840 to square the panel:
# point sizes are absolute while the frame grew ~22%, so unscaled text had quietly shrunk
# relative to everything around it. Raising it further is one edit, not fifteen.
TS <- 1.3

xr <- c(1e-6, 1e4); yr <- c(1e-5, 1)
xs <- 10^seq(log10(xr[1]), log10(xr[2]), length.out = 500)
# The three boundaries as ONE data frame, so they can be a legend instead of a sentence of
# grey text under the plot. They are distinguished by LINETYPE and share a colour, because
# ggplot allows one colour scale per plot and that one is spent on the tiers.
LINE_COL <- "grey30"
# ONE LINE, NOT THREE. Leon, 2026-10-28: the two sloped curves are gone and only the flat
# abundance floor stays. The null y = 1-exp(-x) and the shape test SAT_FLOOR*(1-exp(-x)) were
# drawing the detection geometry across the whole panel, and on a figure whose job is to show
# WHICH TAXA cleared the floor that geometry is backdrop, not content.
#   WHAT THIS COSTS, so it is a choice and not an accident: the saturation floor is one of the
# three tests 05_detect applies, and with its curve gone the figure no longer shows WHY a
# point sitting above the coverage line was still rejected. That information survives only as
# the `sat` column in the side table. If a reviewer asks to see the shape test, it comes back.
# The label is derived from COV_FLOOR so the two cannot drift apart.
LINE_LAB <- c(floor = sprintf("abundance floor (coverage %g%%)", COV_FLOOR * 100))
bounds <- data.frame(x = xs, y = COV_FLOOR, which = "floor")
bounds$which <- factor(bounds$which, levels = names(LINE_LAB))

# SPLIT ONCE. `d[d$run == run, ]` scans all 637,306 rows for every run, which is invisible
# over the six-run sample and quadratic over the cohort: the first full render managed 19
# runs/min against the ~50 the plotting alone costs. Splitting up front trades one pass and
# some memory for a constant-time lookup per run.

# ── the what-if call ─────────────────────────────────────────────────────────────────────────
# MIRRORS detect.py's call_one(): absolute floor, then abundance, then shape, in that order,
# every test POSITIVE so a NA cannot pass by comparing False against everything. It exists so
# the app's sliders can show what a different floor would do WITHOUT re-running the pipeline.
#   The batch renderer does NOT use it: there, `call` comes from detect.py and is authoritative.
# Two implementations of one rule is a known trap, so the app labels its own calls as what-if
# and shows the pipeline's verdict unless a slider has been moved off the default.
recall_rows <- function(s, cov_floor = COV_FLOOR, sat_floor = SAT_FLOOR,
                        min_distinct = 10000) {
  ok <- !is.na(s$distinct) & s$distinct >= min_distinct &
        !is.na(s$coverage) & s$coverage >= cov_floor &
        !is.na(s$saturation) & s$saturation >= sat_floor
  s$call <- ifelse(ok, "yes", "no")
  s$reason <- ifelse(ok, "",
    ifelse(is.na(s$distinct) | s$distinct < min_distinct, "below absolute floor",
    ifelse(s$coverage < cov_floor, "below coverage floor", "reads piled, not spread")))
  s
}

make_panel <- function(s, LABELS = "none", STD_ROWS = 30) {
  if (!nrow(s)) return(NULL)
  run <- s$run[1]
  m <- s[1, ]
  rej <- s[s$call != "yes", ]
  cal <- s[s$call == "yes" & !(s$tier %in% HIDE), ]
  cal$tier <- factor(cal$tier, levels = names(TIER_COLOUR))
  # One mark per row: trophic class for pathogens, the tier itself otherwise. A pathogen whose
  # genus FungalTraits does not carry falls back to the neutral circle rather than being
  # dropped or silently promoted to the strongest class.
  cal$mark <- ifelse(cal$tier == "host", "host",
              ifelse(cal$tier %in% HOLLOW, "resident",
              ifelse(cal$trophic %in% TROPHIC_ORDER, cal$trophic, "necrotroph")))
  cal$mark <- factor(cal$mark, levels = names(MARK_SHAPE))
  cal <- cal[order(-cal$coverage), ]

  # Only reportable pathogens and the host are numbered. The other fungi above the floor are
  # coloured by how they entered the reference and left unlabelled: naming twenty Alternaria
  # and Cladosporium that the Methods does not count as pathogens buries the ones it does.
  cal$named <- cal$tier %in% NAMED
  # NUMBERED IN TIER ORDER, not by coverage alone: the table is read as groups, so the colours
  # cluster and the eye can take "what did PHI-base call on this crop" in one block instead of
  # tracking one hue down an interleaved list. Within a tier it is still coverage, descending.
  cal <- cal[order(cal$tier, -cal$coverage), ]
  # 36 of 2,862 runs clear no floor at all, and that is a RESULT -- the panel still shows
  # the rejected cloud, which is the evidence the floors reject far more than they admit.
  # Assigning into a zero-row frame is what broke the first full render.
  cal$n <- if (nrow(cal)) NA_integer_ else integer(0)
  if (nrow(cal)) cal$n[cal$named] <- seq_len(sum(cal$named))
  # "Zymoseptoria tritici" -> "Z. tritici". Full binomials will not fit beside ten points in
  # the band where everything clusters; the genus initial is unambiguous within one panel.
  cal$abbr <- sub("^([A-Za-z])[A-Za-z]* ", "\\1. ", cal$species)
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
                    tier = factor(names(TIER_COLOUR), levels = names(TIER_COLOUR)),
                    mark = factor("necrotroph", levels = names(MARK_SHAPE)))

  g <- ggplot(s, aes(redundancy, coverage)) +
    geom_line(data = bounds, aes(x, y, linetype = which), inherit.aes = FALSE,
              colour = LINE_COL, linewidth = 0.65) +
    scale_linetype_manual(values = c(floor = "dashed"),
                          labels = LINE_LAB, breaks = names(LINE_LAB), name = NULL) +
    geom_point(data = rej, shape = 21, size = 1.5, stroke = 0.3,
               colour = "grey60", fill = NA, alpha = 0.5) +
    # EVERY tier needs a row somewhere or its legend key renders empty. drop = FALSE keeps the
    # break in the legend but ggplot draws no glyph for a level no layer supplies, which is how
    # "PHI-base: pathogen of another plant" and "congener" came out as labels with no swatch.
    # An invisible point per tier fixes it, and the key stays complete across every run so the
    # legend does not change shape as you page through.
    geom_point(data = key, aes(colour = tier, shape = mark), size = 0,
               fill = NA, alpha = 0, show.legend = TRUE) +
    # SHAPE 16, not 21. A borderless point has to map COLOUR, because shape 21's colour is its
    # border and setting that to NA with stroke 0 makes ggplot drop the glyph entirely -- the
    # first attempt rendered a plot with numbers and no points at all.
    geom_point(data = other, aes(colour = tier, alpha = tier, shape = mark, fill = tier), size = 2.5, stroke = 1.1) +
    geom_point(data = lab, aes(colour = tier, alpha = tier, shape = mark, fill = tier), size = 4.4, stroke = 1.1) +
    # DECLARED is a halo, not a fill, so it composes with any tier colour instead of
    # competing with it. A pathogen the study named carries a ring; everything else does not.
    geom_point(data = cal[cal$declared == "declared", ], shape = 21, size = 7.0,
               stroke = 0.9, colour = INK, fill = NA) +
    (if (LABELS == "none") NULL else
      ggrepel::geom_text_repel(
        data = lab, aes(label = if (LABELS == "name") abbr else as.character(n)),
        size = (if (LABELS == "name") 2.8 else 3.0) * TS, seed = 1,
        inherit.aes = TRUE, colour = INK,
        fontface = if (LABELS == "name") "plain" else "bold", max.overlaps = Inf,
        min.segment.length = 0.1, segment.size = 0.2,
        segment.colour = "grey70", box.padding = 0.3, point.padding = 0.3)) +
    scale_colour_manual(values = TIER_COLOUR, labels = TIER_LABEL,
                        breaks = names(TIER_LABEL), name = NULL, drop = FALSE,
                        aesthetics = "colour") +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none", drop = FALSE) +
    scale_shape_manual(values = MARK_SHAPE, guide = "none", drop = FALSE) +
    scale_fill_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    # override.aes must restore alpha: the invisible key layer above is drawn at alpha = 0 and
    # the legend glyph INHERITS that, so any tier with no real data in a run came out as a
    # label with no swatch. Only visible on runs missing a tier, which is most of them.
    # The legend keys carry the SAME alpha as the points, so the ramp is part of the key
    # rather than being flattened to solid and contradicting the plot. override.aes takes a
    # vector, in the order of `breaks`.
    guides(colour = guide_legend(
             ncol = 2, order = 1,
             override.aes = list(size = 3.6, fill = NA, stroke = 1.1,
                                 shape = ifelse(names(TIER_LABEL) %in% HOLLOW, 1L, 16L),
                                 alpha = unname(TIER_ALPHA[names(TIER_LABEL)]))),
           linetype = guide_legend(ncol = 1, order = 2,
                                   override.aes = list(colour = LINE_COL))) +
    scale_x_log10(limits = c(xs[1], xs[length(xs)]),
                  breaks = 10^seq(-6, 4, by = 2), labels = seq(-6, 4, by = 2)) +
    scale_y_log10(limits = yr, breaks = 10^seq(-5, 0, by = 1),
                  labels = seq(-5, 0, by = 1)) +
    labs(x = expression(log[10]~"Redundancy"), y = expression(log[10]~"Coverage")) +
    theme_classic(base_size = 12 * TS) +
    # SQUARE PANEL. Both axes are decades of the same kind of quantity -- redundancy and
    # coverage are each normalised by the taxon's own minimizer total -- so a decade should be
    # the same distance whichever way it runs. It was not: the canvas height was set by how
    # many rows the side table needed, which left the panel stretched vertically by whatever
    # a run happened to detect. aspect.ratio pins it, and the canvas height below is sized so
    # the square lands without a band of dead space around it.
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.grid       = element_blank(),
          axis.line        = element_line(colour = INK, linewidth = 0.7),
          axis.ticks       = element_line(colour = INK, linewidth = 0.7),
          axis.text        = element_text(colour = INK, size = 13 * TS),
          axis.title       = element_text(colour = INK, size = 15 * TS),
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
  # EVERY tier gets a heading, including the ones that are only counted. Leon, 2026-10-28.
  # The heading carries the swatch and the rows beneath it carry none, so a colour is stated
  # once per block rather than restated on every line; and the resident blocks stop being a
  # bare tally line floating under the named ones and become sections like the rest.
  emit_rows <- function(grp, tr) {
    for (i in seq_len(nrow(grp))) {
      y <<- y - 1
      rows_tb[[length(rows_tb) + 1]] <<- data.frame(
        tier = tr, mark = as.character(grp$mark[i]), y = y,
        txt = if (LABELS == "number") sprintf("%-2d  %s", grp$n[i], grp$species[i])
              else grp$species[i],
        cov = sprintf("%.1f%%", 100 * grp$coverage[i]),
        sat = sprintf("%.2f", grp$saturation[i]),
        declared = grp$declared[i], kind = "row")
    }
  }
  for (tr in NAMED) {
    grp <- cal[cal$tier == tr, ]
    if (!nrow(grp)) next
    y <- y - 1
    rows_tb[[length(rows_tb) + 1]] <- data.frame(
      tier = tr, mark = if (tr == "host") "host" else "necrotroph", y = y,
      txt = TIER_HEAD[[tr]], cov = "", sat = "", declared = "", kind = "head")
    if (tr == "host") {
      emit_rows(grp, tr)
    } else {
      # SUB-HEADED BY TROPHIC CLASS, ordered biotroph -> necrotroph -> opportunist, which is
      # descending strength of claim rather than descending count. The sub-heading carries the
      # SHAPE and the heading above it carries the COLOUR, so each mark is decoded once, at
      # the level it applies to. This is also what pays for dropping the ev column: the table
      # now groups by something a reader can act on instead of by which database said so.
      for (tc in TROPHIC_ORDER) {
        sub_g <- grp[as.character(grp$mark) == tc, ]
        if (!nrow(sub_g)) next
        y <- y - 1
        rows_tb[[length(rows_tb) + 1]] <- data.frame(
          tier = tr, mark = tc, y = y, txt = MARK_LABEL[[tc]],
          cov = "", sat = "", declared = "", kind = "subhead")
        emit_rows(sub_g, tr)
      }
    }
    y <- y - 0.55
  }
  # ONE block for everything that is not a pathogen. The heading carries no swatch, because
  # it spans three colours; the colours sit on the rows under it, which is the only place in
  # the table where the swatch is not on the heading and the only place it would mean
  # anything else.
  hol <- intersect(names(TIER_COLOUR), HOLLOW)
  hol <- hol[sapply(hol, function(tr) any(cal$tier == tr))]
  if (length(hol)) {
    y <- y - 1
    rows_tb[[length(rows_tb) + 1]] <- data.frame(
      tier = hol[1], mark = "resident", y = y, txt = HOLLOW_HEAD, cov = "", sat = "",
      declared = "", kind = "group")
    for (tr in hol) {
      y <- y - 1
      rows_tb[[length(rows_tb) + 1]] <- data.frame(
        tier = tr, mark = "resident", y = y,
        txt = sprintf("%-2d  %s", sum(cal$tier == tr), TIER_ROW[[tr]]),
        cov = "", sat = "", declared = "", kind = "count")
    }
    y <- y - 0.55
  }
  if (!length(rows_tb)) {
    rows_tb[[1]] <- data.frame(
      tier = "no_record", mark = "resident", y = -1,
      txt = "no detection cleared the floors", cov = "", sat = "",
      declared = "", kind = "group")
  }
  tb <- do.call(rbind, rows_tb)
  tb$tier <- factor(tb$tier, levels = names(TIER_COLOUR))
  tb$mark <- factor(tb$mark, levels = names(MARK_SHAPE))

  tg <- ggplot(tb, aes(0, y)) +
    # SWATCH ON THE HEADING, not on the rows -- the reverse of the first draft. The heading
    # is where the colour is defined, so it is where the colour is shown; repeating the dot
    # down every row restated a fact the block header had already made.
    # The block heading carries the COLOUR as a plain filled dot, because the shape belongs
    # to the sub-heading one level below it. Sub-headings and resident rows carry the SHAPE.
    geom_point(data = tb[tb$kind == "head", ], aes(x = 0, colour = tier, alpha = tier),
               shape = 16, size = 3.8) +
    geom_point(data = tb[tb$kind %in% c("subhead", "count"), ],
               aes(x = 0.45, colour = tier, alpha = tier, shape = mark, fill = tier), size = 3.8, stroke = 1.2) +
    # The declared ring keeps a filled centre, so it mirrors the plot exactly: a coloured
    # point inside a ring. It is the only mark left in the body of the table, which suits a
    # figure whose question is which pathogens the study did NOT declare.
    geom_point(data = tb[tb$declared == "declared" & tb$kind == "row", ],
               aes(x = 0.9, colour = tier, shape = mark), size = 3.2, fill = NA,
               show.legend = FALSE) +
    geom_point(data = tb[tb$declared == "declared" & tb$kind == "row", ], aes(x = 0.9),
               shape = 21, size = 6.4, stroke = 0.9, colour = INK, fill = NA) +
    # BOLD ON THE HEADINGS ONLY. The earlier rule was no bold anywhere, because weight was
    # being asked to mark pathogen-versus-host AND named-versus-counted at once and neither
    # reading survived. Weight now marks one thing, section heading against row, and the two
    # readings it used to collide with are carried by hue and by fill.
    geom_text(data = tb[tb$kind %in% c("head", "group"), ], aes(x = 0.35, label = txt),
              hjust = 0, size = 3.1 * TS, colour = INK, fontface = "bold") +
    geom_text(data = tb[tb$kind == "subhead", ], aes(x = 0.8, label = txt), hjust = 0,
              size = 2.85 * TS, colour = "grey35") +
    geom_text(data = tb[tb$kind == "row", ], aes(x = 1.25, label = txt), hjust = 0,
              size = 3.1 * TS, colour = INK) +
    geom_text(data = tb[tb$kind == "count", ], aes(x = 0.8, label = txt), hjust = 0,
              size = 3.1 * TS, colour = "grey45") +
    # THE ev COLUMN IS GONE, Leon 2026-10-28. A/B/C coded which database said so, which is
    # provenance rather than a finding; it belongs in the caption. What replaced it is the
    # trophic sub-heading, which groups by something a reader can act on.
    geom_text(aes(x = 8.8, label = cov), hjust = 1, size = 3.0 * TS, family = "mono",
              colour = INK) +
    geom_text(aes(x = 10.1, label = sat), hjust = 1, size = 3.0 * TS, family = "mono",
              colour = "grey40") +
    annotate("text", x = 8.8, y = 0.1, label = "cov", hjust = 1, size = 2.9 * TS,
             colour = "grey40", fontface = "bold") +
    annotate("text", x = 10.1, y = 0.1, label = "sat", hjust = 1, size = 2.9 * TS,
             colour = "grey40", fontface = "bold") +
    scale_colour_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none", drop = FALSE) +
    scale_shape_manual(values = MARK_SHAPE, guide = "none", drop = FALSE) +
    scale_fill_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    scale_x_continuous(limits = c(-0.5, 10.3)) +
    scale_y_continuous(limits = c(min(-STD_ROWS, min(tb$y)) - 0.8, 0.8)) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(4, 4, 4, 2))

  # WIDTHS sized so the plot column lands square at the height the table dictates, rather
  # than the 1.72 that happened to look balanced when the panel was free to stretch.
  #   The table column is then whatever is left, and it has to be sized in PIXELS against the
  # type size, not chosen by eye: geom_text is absolute while the x scale is data units, so a
  # species name occupies a fraction of the column that moves with TS. At TS = 1.3 the longest
  # name (Cladosporium cladosporioides) ran into the `ev` column at the old 900 px. The column
  # is ~1,180 px here, the same +30% the text gained.
  p <- (g | tg) + patchwork::plot_layout(widths = c(1.661, 1)) +
    patchwork::plot_annotation(
      title = ttl, subtitle = sub,
      theme = theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.title      = element_text(face = "bold", size = 17 * TS),
        plot.subtitle   = element_text(size = 10.5 * TS, colour = "grey25", lineheight = 1.3),
        plot.margin     = margin(12, 14, 8, 12)))

  attr(p, "n_called") <- nrow(cal)
  attr(p, "n_named")  <- nrow(lab)
  attr(p, "n_undeclared") <- n_und
  attr(p, "n_rows") <- max(STD_ROWS, abs(min(tb$y)))
  p
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
  # The eight tier keys are GONE: each table block is now headed by its route name with the
  # swatches beside the rows it applies to, so the colour is decoded where it is used rather
  # than in a key held in the reader's head. What is left is what the table cannot say --
  # the three boundaries, and the ring.
  keep <- setdiff(names(TIER_LABEL), HIDE)
  # The COLOUR keys answer one question, so they are drawn with one mark: filled where the
  # taxon is a pathogen or the host, hollow where it is not. The shapes get their own key
  # below, because conflating them here would imply the two channels are one.
  tl <- data.frame(tier = factor(keep, levels = names(TIER_COLOUR)),
                   lab = unname(TIER_LABEL[keep]),
                   sh = ifelse(keep %in% HOLLOW, 21L, 16L),
                   x = 0, y = -(seq_along(keep) - 1))
  ll <- data.frame(lab = factor(unname(LINE_LAB), levels = unname(LINE_LAB)),
                   x = 7.8, y = -(seq_along(LINE_LAB) - 1))
  dl <- data.frame(lab = "pathogen the study declared", x = 7.8, y = -4.2)
  # The evidence grades, spelled out once. They are a table column now, not a colour, so the
  # legend is the only place the letters are expanded.
  # WHAT KIND OF ORGANISM, the second channel. Drawn in ink rather than in any tier colour,
  # because the shape means the same thing whichever colour it lands in -- which is the whole
  # reason shape and not hue was used for it.
  sk <- c("biotroph", "necrotroph", "opportunist", "host")
  el <- data.frame(sh = unname(MARK_SHAPE[sk]), lab = unname(MARK_LABEL[sk]),
                   x = 7.8, y = -(6.2 + 0.85 * (seq_along(sk) - 1)))

  g <- ggplot() +
    geom_point(data = tl, aes(x, y, colour = tier, alpha = tier, shape = I(sh)),
               fill = NA, size = 4.6, stroke = 1.3) +
    geom_text(data = tl, aes(x + 0.22, y, label = lab), hjust = 0, size = 4.0, colour = INK) +
    geom_segment(data = ll, aes(x = x, xend = x + 0.75, y = y, yend = y, linetype = lab),
                 colour = LINE_COL, linewidth = 0.7) +
    geom_text(data = ll, aes(x + 0.95, y, label = lab), hjust = 0, size = 4.0, colour = INK) +
    # The declared halo. It was missing from the key entirely, which left the one ring on the
    # plot unexplained -- and it is the mark that answers the paper's own question, so it is
    # the last thing that should have gone unlabelled.
    geom_point(data = dl, aes(x + 0.37, y), shape = 16, size = 4.6,
               colour = TIER_COLOUR[["pathogen_crop"]]) +
    geom_point(data = dl, aes(x + 0.37, y), shape = 21, size = 8.4, stroke = 0.9,
               colour = INK, fill = NA) +
    geom_text(data = dl, aes(x + 0.95, y, label = lab), hjust = 0, size = 4.0, colour = INK) +
    annotate("text", x = 7.8, y = -5.4, label = "kind of organism", hjust = 0,
             size = 3.5, colour = "grey45", fontface = "bold") +
    geom_point(data = el, aes(x + 0.37, y, shape = I(sh)), colour = INK, fill = NA,
               size = 4.6, stroke = 1.3) +
    geom_text(data = el, aes(x + 0.95, y, label = lab), hjust = 0, size = 3.8,
              colour = INK) +
    scale_colour_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none", drop = FALSE) +
    scale_shape_manual(values = MARK_SHAPE, guide = "none", drop = FALSE) +
    scale_fill_manual(values = TIER_COLOUR, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = c("dashed"), guide = "none") +
    scale_x_continuous(limits = c(-0.3, 15.0), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-9.6, 0.6), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(14, 14, 14, 14))

  f <- file.path(OUTDIR, "_legend.png")
  ragg::agg_png(f, width = 1900, height = 900, res = 190, background = "white")
  print(g); invisible(dev.off())
  cat("  legend -> _legend.png\n")
}