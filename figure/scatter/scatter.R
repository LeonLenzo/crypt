#!/usr/bin/env Rscript
# figure/scatter/scatter.R — host % vs max pathogen % scatter, faceted by mode
# Run from crypt/: Rscript figure/scatter/scatter.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

RUNS_TSV <- "output/02_filter_runs/data/runs.tsv"
OUT       <- "figure/scatter/scatter"

FLAG_COLS <- c(
  single = "#536d589a",
  coinf  = "#e67e22",
  novel  = "#c02bac"
)
FLAG_LABS <- c(
  single = "Single pathogen",
  coinf  = "Co-infection (known)",
  novel  = "Novel host range"
)

MODE_SHAPES <- c("MAL" = 16, "HAL" = 17)   # circle, triangle

crypt <- read.delim(RUNS_TSV, stringsAsFactors = FALSE) |>
  filter(biosample_representative == "True") |>
  mutate(
    max_pathogen_pct  = fungi_pct + virus_pct + bacteria_pct +
                        oomycete_pct + nematode_pct,
    co_infection_flag = factor(
                         case_when(
                           co_infection_flag == "single"          ~ "single",
                           interaction_status == "novel_host_range" ~ "novel",
                           TRUE                                   ~ "coinf"
                         ),
                         levels = c("single", "coinf", "novel")),
    mode_label        = toupper(mode)
  ) |>
  filter(host_pct > 0, max_pathogen_pct > 0) |>
  arrange(co_infection_flag)   # single drawn first so co-infected sits on top

n_by_flag <- crypt |>
  count(co_infection_flag) |>
  mutate(co_infection_flag = as.character(co_infection_flag))

flag_labs_n <- setNames(
  paste0(FLAG_LABS[n_by_flag$co_infection_flag],
         " (n=", formatC(n_by_flag$n, format = "d", big.mark = ","), ")"),
  n_by_flag$co_infection_flag
)

cat("BioSamples per flag:\n")
print(as.data.frame(n_by_flag))

p <- ggplot(crypt,
            aes(x = host_pct, y = max_pathogen_pct,
                colour = co_infection_flag,
                shape  = mode_label)) +
  geom_point(alpha = 0.3, size = 1.2, stroke = 0) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1%", "1%", "10%", "100%")
  ) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100),
    labels = c("0.1%", "1%", "10%", "100%")
  ) +
  scale_colour_manual(values = FLAG_COLS, labels = flag_labs_n, name = NULL) +
  scale_shape_manual(values = MODE_SHAPES, name = "Mode") +
  guides(
    colour = guide_legend(override.aes = list(alpha = 1, size = 3)),
    shape  = guide_legend(override.aes = list(alpha = 1, size = 3))
  ) +
  labs(
    title    = "Plant host signal vs total pathogen content",
    subtitle = NULL,
    x        = "log₁₀ Viridiplantae kmers (%) ",
    y        = "log₁₀ Pathogen kmers (%) "
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 13),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 17),
    plot.subtitle    = element_text(colour = "grey50", size = 11),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(paste0(OUT, ".pdf"), p, width = 12, height = 6, bg = "white",
       device = cairo_pdf)
ggsave(paste0(OUT, ".png"), p, width = 12, height = 6, bg = "white", dpi = 200)
cat("Written:", paste0(OUT, ".pdf/.png\n"))
