#!/usr/bin/env Rscript
# figure/novel_heatmap/novel_heatmap.R — secondary pathogen × host heatmap for novel_host_range runs
# Run from crypt/: Rscript figure/novel_heatmap/novel_heatmap.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
})

RUNS_TSV    <- "stat/output/filter_runs/data/runs.tsv"
OUT         <- "metadata/output/figures/novel_heatmap/novel_heatmap"
TOP_PATHS   <- 35
TOP_HOSTS   <- 30

crypt <- read.delim(RUNS_TSV, stringsAsFactors = FALSE) |>
  filter(biosample_representative == "True",
         interaction_status       == "novel_host_range",
         stat_pathogens           != "")

cat(sprintf("novel_host_range BioSamples: %d\n", nrow(crypt)))

strip_pct  <- function(x) str_remove_all(x, ":[\\d.]+%")
norm_name  <- function(x) {
  x <- str_trim(x)
  if (str_detect(x, regex("vir(us|oid|inae|ales)", ignore_case = TRUE))) return(x)
  m <- str_match(x, "^(\\S+\\s+\\S+\\s+f\\.\\s+sp\\.\\s+\\S+)")
  if (!is.na(m[,2])) return(m[,2])
  words <- str_split(x, "\\s+")[[1]]
  paste(head(words, 2), collapse = " ")
}

pairs <- crypt |>
  select(host, stat_pathogens, mode) |>
  mutate(
    host               = str_extract(host, "^\\S+\\s+\\S+"),
    stat_pathogens = strip_pct(stat_pathogens)
  ) |>
  separate_rows(stat_pathogens, sep = ";\\s*") |>
  filter(stat_pathogens != "",
         !str_detect(stat_pathogens,
                     regex("environmental|unclassified", ignore_case = TRUE))) |>
  mutate(secondary = sapply(stat_pathogens, norm_name)) |>
  filter(!is.na(host), !is.na(secondary))

top_paths <- pairs |>
  count(secondary, sort = TRUE) |>
  slice_head(n = TOP_PATHS) |>
  pull(secondary)

top_hosts <- pairs |>
  count(host, sort = TRUE) |>
  slice_head(n = TOP_HOSTS) |>
  pull(host)

plot_dat <- pairs |>
  filter(secondary %in% top_paths, host %in% top_hosts) |>
  count(host, secondary) |>
  mutate(
    host      = factor(host,      levels = rev(top_hosts)),
    secondary = factor(secondary, levels = rev(top_paths))
  )

cat(sprintf("Pairs in heatmap: %d  |  hosts: %d  |  pathogens: %d\n",
            nrow(plot_dat), n_distinct(plot_dat$host),
            n_distinct(plot_dat$secondary)))

p <- ggplot(plot_dat, aes(x = secondary, y = host, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n >= 3, n, "")),
            size = 2.8, colour = "white", fontface = "bold") +
  scale_fill_gradient(low = "#fde8d8", high = "#c0392b",
                      name = "BioSamples", na.value = "grey95") +
  scale_x_discrete(position = "bottom") +
  labs(
    title    = "Novel host range detections",
    subtitle = paste0("Secondary pathogen × host species  ·  ",
                      "interaction_status = novel_host_range  ·  ",
                      "one cell per (host, pathogen) pair"),
    x = "Secondary pathogen",
    y = "Host species"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8,
                                    face = "italic"),
    axis.text.y      = element_text(size = 8, face = "italic"),
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey50", size = 9),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(paste0(OUT, ".pdf"), p, width = 18, height = 13, bg = "white",
       device = cairo_pdf)
ggsave(paste0(OUT, ".png"), p, width = 18, height = 13, bg = "white", dpi = 200)
cat("Written:", paste0(OUT, ".pdf/.png\n"))
