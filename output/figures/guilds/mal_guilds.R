#!/usr/bin/env Rscript
# figure/guilds/mal_guilds.R — co-infection guild network (MAL + HAL)
# Run from crypt/: Rscript figure/guilds/mal_guilds.R

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

NODES_TSV <- "figure/guilds/mal_guild_nodes.tsv"
EDGES_TSV <- "figure/guilds/mal_guild_edges.tsv"
OUT_PDF   <- "figure/guilds/mal_guild_network.pdf"
OUT_PNG   <- "figure/guilds/mal_guild_network.png"

MIN_NODE_TOTAL <- 2    # drop nodes never appearing in a kept edge
LABEL_THRESHOLD <- 8   # label nodes with total detections >= this

# ── Colours ───────────────────────────────────────────────────────────────────
KINGDOM_COLS <- c(
  Fungi     = "#d35400",
  Viruses   = "#c0392b",
  Bacteria  = "#2980b9",
  Oomycota  = "#8e44ad",
  Nematoda  = "#27ae60",
  Unknown   = "#95a5a6"
)

# ── Load data ──────────────────────────────────────────────────────────────────
nodes <- read.delim(NODES_TSV, stringsAsFactors = FALSE)
edges <- read.delim(EDGES_TSV, stringsAsFactors = FALSE)

# Keep only nodes that appear in at least one edge
nodes_in_edges <- unique(c(edges$node1, edges$node2))
nodes <- nodes %>% filter(name %in% nodes_in_edges)

# ── Build igraph ───────────────────────────────────────────────────────────────
g <- graph_from_data_frame(
  d        = edges,
  vertices = nodes,
  directed = FALSE
)

# Node attributes
V(g)$kingdom <- nodes$kingdom[match(V(g)$name, nodes$name)]
V(g)$total   <- nodes$total[match(V(g)$name, nodes$name)]
V(g)$label   <- ifelse(nodes$total[match(V(g)$name, nodes$name)] >= LABEL_THRESHOLD,
                       V(g)$name, NA_character_)

# Shorten long virus names for label readability
V(g)$label <- sub("\\s+virus$", " virus", V(g)$label)   # no-op but keeps intent clear
V(g)$label <- gsub("f\\. sp\\. ", "f.sp.", V(g)$label)  # compact f. sp.

# ── Layout ─────────────────────────────────────────────────────────────────────
set.seed(42)
layout <- create_layout(g, layout = "fr")

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggraph(layout) +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour = "#666666",
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.3, 3)) +
  scale_edge_alpha(range = c(0.2, 0.7)) +
  geom_node_point(
    aes(size = total, colour = kingdom),
    alpha = 0.9
  ) +
  geom_node_label(
    aes(label = label),
    size          = 2.2,
    label.padding = unit(0.12, "lines"),
    fill          = "white",
    colour        = "black",
    alpha         = 0.85,
    na.rm         = TRUE,
    repel         = FALSE
  ) +
  scale_colour_manual(values = KINGDOM_COLS, name = "Kingdom") +
  scale_size_continuous(
    name   = "Total detections",
    range  = c(2, 10),
    breaks = c(5, 20, 50, 100, 200)
  ) +
  labs(
    title    = "Cryptic co-infection guilds (MAL + HAL)",
    subtitle = "Nodes = pathogens; edges = co-occurrence in confirmed runs; edge width proportional to shared-run count"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position  = "right",
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold")
  )

ggsave(OUT_PDF, p, width = 12, height = 9)
ggsave(OUT_PNG, p, width = 12, height = 9, dpi = 200)
cat(sprintf("Written: %s, %s\n", OUT_PDF, OUT_PNG))
cat(sprintf("Nodes in graph: %d  |  Edges: %d\n", vcount(g), ecount(g)))
