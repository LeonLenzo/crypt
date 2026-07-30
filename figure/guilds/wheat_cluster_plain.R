#!/usr/bin/env Rscript
# figure/guilds/wheat_cluster_plain.R — plain cereal cluster, no host hulls.
# Uses unfiltered runs (including broad-clade/unresolved hosts).
# Run from crypt/: Rscript figure/guilds/wheat_cluster_plain.R

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

NODES_TSV <- "figure/guilds/field_hc_all_nodes.tsv"
EDGES_TSV <- "figure/guilds/field_hc_all_edges.tsv"
OUT_PDF   <- "figure/guilds/wheat_cluster_plain.pdf"
OUT_PNG   <- "figure/guilds/wheat_cluster_plain.png"

# ── Colour palettes ───────────────────────────────────────────────────────────

ORDER_COLS <- c(
  Pucciniales    = "#c0392b",
  Pleosporales   = "#e67e22",
  Capnodiales    = "#f1c40f",
  Hypocreales    = "#e91e63",
  Glomerellales  = "#9b59b6",
  Magnaporthales = "#3498db",
  Unknown        = "#95a5a6"
)

ORDER_MAP <- c(
  "Alternaria alternata"              = "Pleosporales",
  "Alternaria solani"                 = "Pleosporales",
  "Bipolaris sorokiniana"             = "Pleosporales",
  "Bipolaris zeicola"                 = "Pleosporales",
  "Colletotrichum higginsianum"       = "Glomerellales",
  "Exserohilum turcicum"              = "Pleosporales",
  "Fusarium graminearum"              = "Hypocreales",
  "Leptosphaeria maculans"            = "Pleosporales",
  "Parastagonospora nodorum"          = "Pleosporales",
  "Puccinia graminis f. sp. tritici"  = "Pucciniales",
  "Puccinia striiformis"              = "Pucciniales",
  "Puccinia striiformis f. sp. tritici" = "Pucciniales",
  "Puccinia triticina"                = "Pucciniales",
  "Pyrenophora tritici-repentis"      = "Pleosporales",
  "Pyricularia oryzae"                = "Magnaporthales",
  "Zymoseptoria tritici"              = "Capnodiales"
)

EXCLUDE_NODES <- c("Listeria monocytogenes", "Plasmopara viticola")

# ── Layout helper ─────────────────────────────────────────────────────────────

stretch_center <- function(layout, power = 0.60) {
  xy     <- as.matrix(layout[, c("x", "y")])
  center <- colMeans(xy)
  delta  <- sweep(xy, 2, center)
  dist   <- sqrt(rowSums(delta^2))
  dist[dist < 1e-9] <- 1e-9
  scale  <- dist^(power - 1)
  xy_new <- sweep(delta * scale, 2, center, "+")
  layout$x <- xy_new[, 1]
  layout$y <- xy_new[, 2]
  layout
}

# ── Load & build graph ────────────────────────────────────────────────────────

nodes <- read.delim(NODES_TSV, stringsAsFactors = FALSE)
edges <- read.delim(EDGES_TSV, stringsAsFactors = FALSE)

nodes_in_edges <- unique(c(edges$node1, edges$node2))
nodes <- nodes |> filter(name %in% nodes_in_edges, !name %in% EXCLUDE_NODES)
edges <- edges |> filter(!node1 %in% EXCLUDE_NODES, !node2 %in% EXCLUDE_NODES)

nodes$order <- ORDER_MAP[nodes$name]
nodes$order[is.na(nodes$order)] <- "Unknown"

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# ── Largest connected component ───────────────────────────────────────────────

comps   <- components(g)
biggest <- which.max(comps$csize)
g_wheat <- induced_subgraph(g, which(comps$membership == biggest))

cat(sprintf("Largest component: %d nodes, %d edges\n", vcount(g_wheat), ecount(g_wheat)))

shorten <- function(x) gsub("f\\. sp\\. ", "f.sp.", x)
V(g_wheat)$label <- shorten(V(g_wheat)$name)
V(g_wheat)$order <- nodes$order[match(V(g_wheat)$name, nodes$name)]

# ── Layout ─────────────────────────────────────────────────────────────────────

set.seed(42)
layout <- create_layout(g_wheat, layout = "fr")
layout <- stretch_center(layout, power = 0.60)

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggraph(layout) +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour      = "#555555",
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.5, 5),  guide = "none") +
  scale_edge_alpha(range = c(0.55, 0.90), guide = "none") +
  geom_node_point(
    aes(size = total, fill = order),
    shape = 21, colour = "white", stroke = 0.4, alpha = 0.95
  ) +
  geom_node_label(
    aes(label = label, colour = order),
    size          = 3.0,
    label.padding = unit(0.15, "lines"),
    fill          = "white",
    alpha         = 0.88,
    repel         = TRUE,
    max.overlaps  = 40,
    show.legend   = FALSE
  ) +
  scale_fill_manual(values = ORDER_COLS, name = "Order",
                    guide = guide_legend(order = 1,
                      override.aes = list(size = 4, shape = 21,
                                          stroke = 0.4))) +
  scale_colour_manual(values = ORDER_COLS, guide = "none") +
  scale_size_continuous(
    name   = "Detections\n(field HC runs)",
    range  = c(3, 14),
    breaks = c(5, 15, 30, 60, 80),
    guide  = guide_legend(order = 2)
  ) +
  labs(
    title    = "Co-infection guild of cereal foliar pathogens in field RNA-seq data",
    subtitle = paste0(
      "Each edge connects pathogens co-detected in the same field biosample",
      " (NCBI SRA); edge width proportional to co-occurrence frequency;\n",
      "node size = detection count; node colour = fungal order",
      " (all host resolution levels included)"
    )
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 13, face = "bold"),
    plot.subtitle   = element_text(size = 8.5, colour = "grey40"),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10, face = "bold")
  )

ggsave(OUT_PDF, p, width = 12, height = 9)
ggsave(OUT_PNG, p, width = 12, height = 9, dpi = 250)
cat(sprintf("Written: %s\n", OUT_PNG))
