#!/usr/bin/env Rscript
# crypt_host_tree.R — plant host tree for MAL confirmed co-infection runs
# Run from crypt/ : Rscript figure/host_tree/crypt_host_tree.R

suppressPackageStartupMessages({
  library(ggtree); library(ape); library(ggplot2); library(dplyr); library(tibble)
  library(cowplot)
})

NWK      <- "metadata/output/figures/host_tree/crypt_host_tree.nwk"
META     <- "metadata/output/figures/host_tree/crypt_host_tree_meta.tsv"
OUT      <- "metadata/output/figures/host_tree/crypt_host_tree"
MIN_RUNS <- 1      # minimum confirmed runs to include a tip
BAR_UNIT  <- 7     # x-units per log10 run unit
BAR_LW    <- 2.5   # linewidth for each bar (mm)
Y_OFF     <- 0.1752  # y-offset: ≈ half of BAR_LW in data units → bars touch

COL_MULTI  <- "#e67e22"   # orange — co-infections (multi_species / multi_kingdom)
COL_SINGLE <- "#27ae60"   # green  — single pathogen detected

# ── Load + filter ─────────────────────────────────────────────────────────────
tree <- read.tree(NWK)
meta <- read.delim(META, stringsAsFactors = FALSE) %>%
  mutate(n_single    = as.integer(n_single),
         n_multi     = as.integer(n_multi),
         n_confirmed = as.integer(n_confirmed),
         family      = ifelse(family == "", NA_character_, family))

keep <- meta %>%
  filter(n_confirmed >= MIN_RUNS) %>%
  pull(label)
tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
n_tips <- Ntip(tree)
cat(sprintf("Tips after filter (>=%d confirmed): %d\n", MIN_RUNS, n_tips))

# ── Base tree ─────────────────────────────────────────────────────────────────
p <- ggtree(tree, layout = "fan", open.angle = 20,
            branch.length = "none",
            size = 0.25, color = "grey45") %<+% meta

# ── Coordinates ───────────────────────────────────────────────────────────────
tree_data <- p$data
max_tip   <- max(tree_data$x[tree_data$isTip], na.rm = TRUE)

BAR_RING  <- max_tip + 1
EXPONENTS <- 0:3
GRID_X    <- BAR_RING + EXPONENTS * BAR_UNIT
GRID_LABS <- as.character(EXPONENTS)
BAND_X    <- max(GRID_X) + 3
NUM_X     <- BAND_X + 2.5
MAX_X     <- NUM_X  + 2

tip_dat <- tree_data %>%
  filter(isTip, !is.na(n_confirmed), n_confirmed >= 1) %>%
  mutate(
    bar_multi_end  = BAR_RING + log10(pmax(n_multi,   1)) * BAR_UNIT,
    bar_single_end = BAR_RING + log10(pmax(n_single,  1)) * BAR_UNIT
  )

# ── Layers ────────────────────────────────────────────────────────────────────
for (gx in GRID_X) {
  p <- p + geom_vline(xintercept = gx, color = "grey85", linewidth = 0.5)
}

p <- p + geom_segment(
  data      = tip_dat,
  aes(x = x, xend = BAR_RING, y = y, yend = y),
  color     = "grey80", linewidth = 0.35, na.rm = TRUE
)

p <- p + geom_segment(
  data      = tip_dat %>% filter(n_multi > 0),
  aes(x = BAR_RING, xend = bar_multi_end, y = y - Y_OFF, yend = y - Y_OFF),
  color     = COL_MULTI, linewidth = BAR_LW, na.rm = TRUE
)

p <- p + geom_segment(
  data      = tip_dat %>% filter(n_single > 0),
  aes(x = BAR_RING, xend = bar_single_end, y = y + Y_OFF, yend = y + Y_OFF),
  color     = COL_SINGLE, linewidth = BAR_LW, na.rm = TRUE
)

p <- p + annotate(
  "text", x = GRID_X, y = -8, label = GRID_LABS,
  size = 9, color = "grey55", hjust = 0.5, vjust = 0.5
)

# ── Rotate ────────────────────────────────────────────────────────────────────
p_rot <- rotate_tree(p, 10)
p_rot <- p_rot + scale_y_continuous(limits = c(-15, NA))

# ── Family bands + letter labels ──────────────────────────────────────────────
rot_data <- p_rot$data

fam_ranges <- rot_data %>%
  filter(isTip) %>%
  select(label, y) %>%
  left_join(meta %>% select(label, family), by = "label") %>%
  filter(!is.na(family)) %>%
  group_by(family) %>%
  summarise(
    y_min      = min(y),
    y_max      = max(y),
    y_mid      = (min(y) + max(y)) / 2,
    n_fam_tips = n(),
    .groups    = "drop"
  ) %>%
  filter(n_fam_tips >= 2) %>%
  arrange(y_mid) %>%
  mutate(fam_lbl = LETTERS[row_number()])

p_rot <- p_rot +
  geom_segment(
    data      = fam_ranges,
    aes(x = BAND_X, xend = BAND_X, y = y_min, yend = y_max),
    color = "#4e4d4d", linewidth = 5, na.rm = TRUE
  ) +
  geom_text(
    data  = fam_ranges,
    aes(x = NUM_X, y = y_mid, label = fam_lbl),
    angle = 0, hjust = 0.5, color = "grey15",
    size = 10, fontface = "plain", na.rm = TRUE
  )

# ── Scales + theme ────────────────────────────────────────────────────────────
total_confirmed <- sum(meta$n_confirmed, na.rm = TRUE)

p_rot <- p_rot +
  xlim(0, MAX_X) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 36, face = "bold",  hjust = 0.5,
                                   color = "grey15"),
    plot.subtitle   = element_text(size = 24, face = "plain", hjust = 0.5,
                                   color = "grey55", margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(20, 5, 20, 20)
  ) +
  labs(
    title    = "Co-infection Candidates by Plant Host",
    
  )

# ── Family legend ─────────────────────────────────────────────────────────────
n_fam   <- nrow(fam_ranges)
ROW_SEP <- 0.05

leg_dat <- fam_ranges %>%
  select(fam_lbl, family) %>%
  mutate(y = rev(seq_len(n_fam)))

p_leg <- ggplot(leg_dat, aes(y = y)) +
  geom_text(aes(x = 0,   label = fam_lbl),
            hjust = 0, size = 10, fontface = "plain", color = "grey15") +
  geom_text(aes(x = 0.5, label = family),
            hjust = 0, size = 8,  fontface = "plain", color = "grey15") +
  scale_x_continuous(limits = c(-0.1, 4)) +
  scale_y_continuous(limits = c(0.5, n_fam + 0.5)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin     = margin(0, 0, 0, 0))

leg_h   <- n_fam * ROW_SEP
leg_bot <- (1 - leg_h) / 2

# ── Combine and save ──────────────────────────────────────────────────────────
combined <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_plot(p_rot, x = 0,    y = 0,       width = 0.72, height = 1) +
  draw_plot(p_leg, x = 0.65, y = leg_bot, width = 0.27, height = leg_h)

ggsave(paste0(OUT, ".pdf"), combined, width = 24, height = 12,
       bg = "white", device = cairo_pdf)
ggsave(paste0(OUT, ".png"), combined, width = 24, height = 12,
       bg = "white", dpi = 200)
cat("Written:", paste0(OUT, ".pdf/.png\n"))
