#!/usr/bin/env Rscript
# study_design_figures.R
# Two figures exploring co-infection rates by LLM-classified study type.
# Run from crypt/: Rscript figure/study_design/study_design_figures.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

OUTDIR <- "output/figures/study_design"

# ── Shared theme ──────────────────────────────────────────────────────────────

base_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

SET_COLS <- c("Field" = "#27ae60", "Lab" = "#8e44ad", "Unclear" = "#b0b0b0")

TREAT_ORDER  <- c("single", "host_study", "abiotic_stress", "surveillance",
                   "coinf_experiment", "combined_stress")
TREAT_LABELS <- c("single\npathogen", "host\nstudy", "abiotic\nstress",
                   "surveillance", "coinf\nexperiment", "combined\nstress")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — Co-infection rate by treatment × setting (all vs high-confidence)
# ══════════════════════════════════════════════════════════════════════════════

d1 <- read.delim("output/figures/study_design/treat_setting_rates.tsv", stringsAsFactors = FALSE)

# Exclude tiny cells and unclear treatment
d1 <- d1 %>%
  filter(treatment %in% TREAT_ORDER, setting %in% c("field", "lab", "unclear"),
         n_bp >= 5) %>%
  mutate(
    treatment = factor(treatment, levels = TREAT_ORDER, labels = TREAT_LABELS),
    setting   = factor(setting,   levels = c("field", "lab", "unclear"),
                                  labels = c("Field", "Lab", "Unclear"))
  )

# Long format for all vs hc
d1_long <- d1 %>%
  select(treatment, setting, n_bp, coinf_rate, coinf_rate_hc) %>%
  pivot_longer(cols = c(coinf_rate, coinf_rate_hc),
               names_to = "filter", values_to = "rate") %>%
  mutate(filter = factor(filter,
                         levels = c("coinf_rate", "coinf_rate_hc"),
                         labels = c("All co-infections",
                                    "High confidence (diff-genus only)")))

p1 <- ggplot(d1_long, aes(x = treatment, y = rate, fill = setting)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = NA) +
  facet_wrap(~ filter, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.06))) +
  scale_fill_manual(values = SET_COLS, name = "Setting") +
  labs(
    title    = "Co-infection rate by study type and setting",
    subtitle = "LLM classification; cells with < 5 BioProjects excluded",
    x        = NULL,
    y        = "Co-infection rate (biosample-rep runs)"
  ) +
  base_theme

ggsave(file.path(OUTDIR, "treat_setting_bar.pdf"), p1, width = 8.5, height = 7)
ggsave(file.path(OUTDIR, "treat_setting_bar.png"), p1, width = 8.5, height = 7, dpi = 200)
cat("Wrote treat_setting_bar.pdf/png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — Single-treatment BPs: field vs lab, all vs high-confidence
# Two panels:
#   A  proportion of BPs with ANY co-infection (all vs hc)
#   B  distribution of coinf_rate among BPs with at least one coinf run
# ══════════════════════════════════════════════════════════════════════════════

d2 <- read.delim("output/figures/study_design/single_bp_rates.tsv", stringsAsFactors = FALSE) %>%
  filter(n_runs >= 2) %>%
  mutate(setting = factor(setting, levels = c("field", "lab"),
                          labels = c("Field", "Lab")))

# Panel A: proportion with any co-infection ──────────────────────────────────

prop_df <- d2 %>%
  group_by(setting) %>%
  summarise(
    n          = n(),
    any_all    = mean(coinf_rate    > 0),
    any_hc     = mean(coinf_rate_hc > 0),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(any_all, any_hc),
               names_to = "filter", values_to = "prop") %>%
  mutate(
    filter = factor(filter,
                    levels = c("any_all", "any_hc"),
                    labels = c("All", "High confidence")),
    label  = percent(prop, accuracy = 1)
  )

pA <- ggplot(prop_df, aes(x = setting, y = prop, fill = setting, alpha = filter)) +
  geom_col(position = position_dodge(0.7), width = 0.6, colour = NA) +
  geom_text(aes(label = label), position = position_dodge(0.7),
            vjust = -0.4, size = 3, fontface = "bold") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 0.75),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = SET_COLS, guide = "none") +
  scale_alpha_manual(values = c("All" = 1, "High confidence" = 0.55),
                     name = NULL) +
  labs(title = "A — % BioProjects with any co-infection",
       x = NULL, y = "Proportion of BioProjects") +
  base_theme +
  theme(legend.position = "bottom")

# Panel B: distribution among BPs with any co-infection ──────────────────────

# Long format: one row per BP per filter type, only those with rate > 0 in either
d2_any <- d2 %>%
  filter(coinf_rate > 0 | coinf_rate_hc > 0) %>%
  select(bp, setting, coinf_rate, coinf_rate_hc) %>%
  pivot_longer(cols = c(coinf_rate, coinf_rate_hc),
               names_to = "filter", values_to = "rate") %>%
  filter(rate > 0) %>%
  mutate(filter = factor(filter,
                         levels = c("coinf_rate", "coinf_rate_hc"),
                         labels = c("All", "High confidence")))

# Median labels
med_df <- d2_any %>%
  group_by(setting, filter) %>%
  summarise(med = median(rate), n = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", n, "\nmed=", percent(med, accuracy = 1)))

pB <- ggplot(d2_any, aes(x = interaction(filter, setting), y = rate,
                          colour = setting, fill = setting, alpha = filter)) +
  geom_violin(trim = TRUE, linewidth = 0.3, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.7, outlier.alpha = 0.4,
               fill = "white", linewidth = 0.35, colour = "grey40") +
  geom_text(data = med_df,
            aes(x = interaction(filter, setting), y = 1.03, label = label),
            size = 2.6, colour = "grey30", vjust = 0, inherit.aes = FALSE) +
  scale_x_discrete(labels = rep(c("All", "High conf."), 2)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.18),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_manual(values   = SET_COLS, guide = "none") +
  scale_colour_manual(values = SET_COLS, guide = "none") +
  scale_alpha_manual(values  = c("All" = 0.65, "High confidence" = 0.45), guide = "none") +
  facet_wrap(~ setting, scales = "free_x") +
  labs(title = "B — Co-infection rate distribution (BPs with ≥ 1 coinf run)",
       x = NULL, y = "Per-BioProject co-infection rate") +
  base_theme +
  theme(strip.text = element_text(colour = c("#27ae60", "#8e44ad")))

p2 <- pA / pB +
  plot_annotation(
    title    = "Single-pathogen BioProjects: field vs lab co-infection profile",
    subtitle = "BioProjects with ≥ 2 biosample-representative runs; LLM treatment = 'single'"
  )

ggsave(file.path(OUTDIR, "single_field_lab.pdf"), p2, width = 8, height = 9)
ggsave(file.path(OUTDIR, "single_field_lab.png"), p2, width = 8, height = 9, dpi = 200)
cat("Wrote single_field_lab.pdf/png\n")
