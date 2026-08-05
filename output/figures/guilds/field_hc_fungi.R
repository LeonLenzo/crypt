#!/usr/bin/env Rscript
# figure/guilds/field_hc_fungi.R — fungi-only co-infection guild network,
# field samples, high confidence (diff-genus secondaries).
# Run from crypt/: Rscript figure/guilds/field_hc_fungi.R

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(ggrepel)

NODES_TSV <- "figure/guilds/field_hc_nodes.tsv"
EDGES_TSV <- "figure/guilds/field_hc_edges.tsv"
OUT_PDF   <- "figure/guilds/field_hc_fungi_network.pdf"
OUT_PNG   <- "figure/guilds/field_hc_fungi_network.png"

LABEL_THRESHOLD <- 3

FUNGI_COL <- "#d35400"

# ── Load + filter to fungi ─────────────────────────────────────────────────────

nodes <- read.delim(NODES_TSV, stringsAsFactors = FALSE) %>%
  filter(kingdom == "Fungi")

fungi_names <- nodes$name

edges <- read.delim(EDGES_TSV, stringsAsFactors = FALSE) %>%
  filter(node1 %in% fungi_names & node2 %in% fungi_names)

nodes <- nodes %>% filter(name %in% unique(c(edges$node1, edges$node2)))

cat(sprintf("Fungi nodes: %d  |  Edges: %d\n", nrow(nodes), nrow(edges)))

# ── Build igraph ───────────────────────────────────────────────────────────────

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
V(g)$total <- nodes$total[match(V(g)$name, nodes$name)]

shorten <- function(x) gsub("f\\. sp\\. ", "f.sp.", x)
V(g)$label <- ifelse(
  nodes$total[match(V(g)$name, nodes$name)] >= LABEL_THRESHOLD,
  shorten(V(g)$name),
  NA_character_
)

# ── Layout + plot ──────────────────────────────────────────────────────────────

set.seed(42)
layout <- create_layout(g, layout = "fr")

p <- ggraph(layout) +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour = "#888888",
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.4, 5)) +
  scale_edge_alpha(range = c(0.2, 0.8)) +
  geom_node_point(
    aes(size = total),
    colour = FUNGI_COL,
    alpha  = 0.90
  ) +
  geom_node_label(
    aes(label = label),
    colour        = FUNGI_COL,
    size          = 3.2,
    label.padding = unit(0.13, "lines"),
    fill          = "white",
    alpha         = 0.90,
    na.rm         = TRUE,
    repel         = TRUE,
    max.overlaps  = 40,
    show.legend   = FALSE
  ) +
  scale_size_continuous(
    name   = "Detections\n(field HC runs)",
    range  = c(2, 14),
    breaks = c(5, 15, 30, 60, 80)
  ) +
  labs(
    title    = "Fungal co-infection guilds — field samples (high confidence)",
    subtitle = paste0(
      "Nodes = fungal pathogens; edges = co-occurrence in same biosample-representative run\n",
      "Filtered to: field (LLM), fungi only, diff-genus secondaries, edge weight >= 1"
    )
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 8,  colour = "grey40"),
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9,  face = "bold")
  )

ggsave(OUT_PDF, p, width = 13, height = 10)
ggsave(OUT_PNG, p, width = 13, height = 10, dpi = 250)
cat(sprintf("Written: %s\n", OUT_PNG))
