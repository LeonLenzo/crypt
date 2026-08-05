#!/usr/bin/env Rscript
# figure/guilds/field_hc_guilds.R — co-infection guild network for high-confidence
# field samples (same_genus_secondary=False, llm_study_setting=field).
# Run from crypt/: Rscript figure/guilds/field_hc_guilds.R

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(ggrepel)

NODES_TSV <- "output/figures/guilds/field_hc_nodes.tsv"
EDGES_TSV <- "output/figures/guilds/field_hc_edges.tsv"
OUT_PDF   <- "output/figures/guilds/field_hc_network.pdf"
OUT_PNG   <- "output/figures/guilds/field_hc_network.png"

LABEL_THRESHOLD <- 4   # label nodes with total detections >= this

KINGDOM_COLS <- c(
  Fungi    = "#d35400",
  Viruses  = "#c0392b",
  Bacteria = "#2980b9",
  Oomycota = "#8e44ad",
  Nematoda = "#27ae60",
  Unknown  = "#95a5a6"
)

# ── Load ───────────────────────────────────────────────────────────────────────

nodes <- read.delim(NODES_TSV, stringsAsFactors = FALSE)
edges <- read.delim(EDGES_TSV, stringsAsFactors = FALSE)

nodes_in_edges <- unique(c(edges$node1, edges$node2))
nodes <- nodes %>% filter(name %in% nodes_in_edges)

# ── Build igraph ───────────────────────────────────────────────────────────────

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

V(g)$kingdom <- nodes$kingdom[match(V(g)$name, nodes$name)]
V(g)$total   <- nodes$total[match(V(g)$name, nodes$name)]

# Shorten long names for legibility
shorten <- function(x) {
  x <- gsub("f\\. sp\\. ", "f.sp.", x)
  x <- gsub(" virus$", " v.", x)
  x
}
V(g)$label <- ifelse(
  nodes$total[match(V(g)$name, nodes$name)] >= LABEL_THRESHOLD,
  shorten(V(g)$name),
  NA_character_
)

# ── Layout ─────────────────────────────────────────────────────────────────────

set.seed(42)
layout <- create_layout(g, layout = "fr")

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggraph(layout) +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour = "#555555",
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.4, 4)) +
  scale_edge_alpha(range = c(0.25, 0.75)) +
  geom_node_point(
    aes(size = total, colour = kingdom),
    alpha = 0.92
  ) +
  geom_node_label(
    aes(label = label, colour = kingdom),
    size          = 3.0,
    label.padding = unit(0.13, "lines"),
    fill          = "white",
    alpha         = 0.90,
    na.rm         = TRUE,
    repel         = TRUE,
    max.overlaps  = 30,
    show.legend   = FALSE
  ) +
  scale_colour_manual(values = KINGDOM_COLS, name = "Kingdom") +
  scale_size_continuous(
    name   = "Detections\n(field HC runs)",
    range  = c(2, 12),
    breaks = c(5, 15, 30, 60, 80)
  ) +
  labs(
    title    = "Cryptic co-infection guilds - field samples (high confidence)",
    subtitle = paste0(
      "Nodes = pathogens; edges = co-occurrence in same biosample-representative run\n",
      "Filtered to: field (LLM), diff-genus secondaries, edge weight >= 2"
    )
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 8, colour = "grey40"),
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9, face = "bold")
  )

ggsave(OUT_PDF, p, width = 14, height = 11)
ggsave(OUT_PNG, p, width = 14, height = 11, dpi = 250)
cat(sprintf("Written: %s\n", OUT_PNG))
cat(sprintf("Nodes: %d  |  Edges: %d\n", vcount(g), ecount(g)))
