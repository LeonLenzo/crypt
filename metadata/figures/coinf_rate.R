#!/usr/bin/env Rscript
# figure/coinf_rate/coinf_rate.R — confirmed BioSamples per host: single vs co-infected (log scale)
# Run from crypt/: Rscript figure/coinf_rate/coinf_rate.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(stringr)
})

RUNS_TSV  <- "stat/output/filter_runs/data/runs.tsv"
OUT        <- "metadata/output/figures/coinf_rate/coinf_rate"
MIN_BS     <- 10
N_HOSTS    <- 20
BAR_W      <- 0.32
GAP        <- 0.06
SPACING    <- 1.8    # y_pos multiplier — increase to widen host gaps

COL_SINGLE <- "#27ae60"
COL_COINF  <- "#e67e22"

host_meta <- read.delim("metadata/output/figures/host_tree/crypt_host_tree_meta.tsv",
                        stringsAsFactors = FALSE) |>
  select(species, family)

crypt <- read.delim(RUNS_TSV, stringsAsFactors = FALSE) |>
  filter(biosample_representative == "True") |>
  mutate(host_sp  = str_extract(host, "^\\S+\\s+\\S+"),
         is_coinf = co_infection_flag != "single") |>
  inner_join(host_meta, by = c("host_sp" = "species"))

count_dat <- crypt |>
  filter(!is.na(host_sp)) |>
  group_by(host_sp) |>
  summarise(n_total  = n(),
            n_coinf  = sum(is_coinf),
            n_single = sum(!is_coinf),
            .groups  = "drop") |>
  filter(n_total >= MIN_BS) |>
  arrange(desc(n_coinf)) |>
  slice_head(n = N_HOSTS) |>
  mutate(y_pos = row_number() * SPACING)

cat(sprintf("Hosts shown: %d  (min %d BioSamples)\n", nrow(count_dat), MIN_BS))

rect_dat <- bind_rows(
  count_dat |>
    transmute(y_pos,
              xmin = y_pos - BAR_W, xmax = y_pos - GAP,
              ymin = 0.5, ymax = n_single + 0.5,
              n    = n_single,
              type = "Single infection"),
  count_dat |>
    filter(n_coinf > 0) |>
    transmute(y_pos,
              xmin = y_pos + GAP, xmax = y_pos + BAR_W,
              ymin = 0.5, ymax = n_coinf + 0.5,
              n    = n_coinf,
              type = "Co-infection")
) |>
  mutate(type = factor(type, levels = c("Single infection", "Co-infection")))

TYPE_COLS <- c("Single infection" = COL_SINGLE, "Co-infection" = COL_COINF)

p <- ggplot(rect_dat) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax, fill = type)) +
  geom_text(aes(x     = (xmin + xmax) / 2,
                y     = ymax,
                label = formatC(n, format = "d", big.mark = ","),
                colour = type),
            hjust = -0.1, vjust = 0.5, angle = 0,
            size = 4.5, fontface = "bold") +
  scale_fill_manual(values = TYPE_COLS, name = NULL) +
  scale_colour_manual(values = TYPE_COLS, guide = "none") +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("1", "10", "100", "1k", "10k"),
    expand = expansion(mult = c(0, 0.3))
  ) +
  scale_x_continuous(
    breaks = count_dat$y_pos,
    labels = count_dat$host_sp,
    expand = expansion(add = c(0.8, 0.8))
  ) +
  coord_flip() +
  labs(
    title    = "Confirmed BioSamples per host species",
    subtitle = paste0("Top ", N_HOSTS, " hosts by total confirmed BioSamples  ·  ",
                      "both modes  ·  log scale"),
    x = NULL,
    y = "BioSamples (log₁₀)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y        = element_text(face = "italic", size = 12),
    axis.text.x        = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 13),
    plot.title         = element_text(face = "bold", size = 16),
    plot.subtitle      = element_text(colour = "grey50", size = 10),
    plot.background    = element_rect(fill = "white", colour = NA),
    plot.margin        = margin(10, 40, 10, 10)
  )

ggsave(paste0(OUT, ".pdf"), p, width = 14, height = 14, bg = "white",
       device = cairo_pdf)
ggsave(paste0(OUT, ".png"), p, width = 14, height = 14, bg = "white", dpi = 200)
cat("Written:", paste0(OUT, ".pdf/.png\n"))
