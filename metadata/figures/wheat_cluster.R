#!/usr/bin/env Rscript
# figure/guilds/wheat_cluster.R — render the wheat foliar pathogen sub-network.
# Extracts the largest connected component from the field HC guild graph.
# Run from crypt/: Rscript figure/guilds/wheat_cluster.R

library(igraph)
library(ggraph)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(dplyr)
library(graphlayouts)

NODES_TSV      <- "metadata/output/figures/guilds/field_hc_nodes.tsv"
EDGES_TSV      <- "metadata/output/figures/guilds/field_hc_edges.tsv"
NODE_HOSTS_TSV <- "metadata/output/figures/guilds/field_hc_node_hosts.tsv"
OUT_PDF        <- "metadata/output/figures/guilds/wheat_cluster.pdf"
OUT_PNG        <- "metadata/output/figures/guilds/wheat_cluster.png"

# Minimum co-infection run count to include a (pathogen, host) pair in a hull.
MIN_HOST_RUNS  <- 1
# Minimum number of distinct cluster nodes a host must span to draw a hull.
MIN_HOST_NODES <- 3
# Whitelist of biologically credible hosts for this cereal network.
# Quercus, Selaginella etc. that pass the numeric threshold are STAT noise.
HULL_WHITELIST <- c(
  "Triticum sp.", "Hordeum vulgare", "Zea mays", "Sorghum bicolor"
)

# ── Colour palettes ───────────────────────────────────────────────────────────

# Fungal order colours (NCBI taxonomy)
ORDER_COLS <- c(
  Pucciniales    = "#c0392b",
  Pleosporales   = "#e67e22",
  Capnodiales    = "#f1c40f",
  Hypocreales    = "#e91e63",
  Glomerellales  = "#9b59b6",
  Magnaporthales = "#3498db",
  Helotiales     = "#27ae60",
  Eurotiales     = "#16a085",
  Agaricales     = "#795548",
  Unknown        = "#95a5a6"
)

ORDER_MAP <- c(
  # Pucciniales (rusts)
  "Puccinia graminis f. sp. tritici"    = "Pucciniales",
  "Puccinia striiformis"                = "Pucciniales",
  "Puccinia striiformis f. sp. tritici" = "Pucciniales",
  "Puccinia triticina"                  = "Pucciniales",
  # Pleosporales
  "Alternaria alternata"                = "Pleosporales",
  "Alternaria solani"                   = "Pleosporales",
  "Bipolaris maydis"                    = "Pleosporales",
  "Bipolaris sorokiniana"               = "Pleosporales",
  "Bipolaris zeicola"                   = "Pleosporales",
  "Cercospora beticola"                 = "Pleosporales",
  "Cercospora kikuchii"                 = "Pleosporales",
  "Exserohilum turcica"                 = "Pleosporales",
  "Exserohilum turcicum"                = "Pleosporales",
  "Leptosphaeria maculans"              = "Pleosporales",
  "Parastagonospora nodorum"            = "Pleosporales",
  "Pyrenophora tritici-repentis"        = "Pleosporales",
  "Setosphaeria turcica"                = "Pleosporales",
  "Stemphylium lycopersici"             = "Pleosporales",
  # Capnodiales
  "Zymoseptoria tritici"                = "Capnodiales",
  # Hypocreales
  "Fusarium graminearum"                = "Hypocreales",
  "Fusarium odoratissimum"              = "Hypocreales",
  "Fusarium oxysporum"                  = "Hypocreales",
  "Fusarium solani"                     = "Hypocreales",
  "Fusarium verticillioides"            = "Hypocreales",
  # Glomerellales
  "Colletotrichum fructicola"           = "Glomerellales",
  "Colletotrichum higginsianum"         = "Glomerellales",
  "Colletotrichum siamense"             = "Glomerellales",
  # Magnaporthales
  "Pyricularia oryzae"                  = "Magnaporthales",
  # Helotiales
  "Botrytis cinerea"                    = "Helotiales",
  "Sclerotinia sclerotiorum"            = "Helotiales",
  # Eurotiales
  "Penicillium expansum"                = "Eurotiales",
  # Agaricales (Basidiomycota)
  "Moniliophthora perniciosa"           = "Agaricales"
)

# Host hull colours — one per host species.  Any unlisted host gets grey.
HOST_PALETTE <- c(
  "Triticum sp."    = "#ffe8a0",
  "Hordeum vulgare" = "#d4b896",
  "Zea mays"        = "#90d890",
  "Sorghum bicolor" = "#e8b090"
)

EXCLUDE_NODES <- c("Listeria monocytogenes", "Plasmopara viticola")

# ── Centre-stretching layout helper ───────────────────────────────────────────
# Applies a power < 1 to node distances from centroid: compresses the outer
# ring and opens up the crowded hub. power=0.6 gives a moderate stretch.

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

fungal_names <- nodes$name[nodes$kingdom == "Fungi"]
nodes <- nodes |> filter(kingdom == "Fungi", !name %in% EXCLUDE_NODES)
edges <- edges |> filter(node1 %in% fungal_names, node2 %in% fungal_names,
                         !node1 %in% EXCLUDE_NODES, !node2 %in% EXCLUDE_NODES)

nodes$order <- ORDER_MAP[nodes$name]
nodes$order[is.na(nodes$order)] <- "Unknown"

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# ── Largest connected component ───────────────────────────────────────────────

comps   <- components(g)
biggest <- which.max(comps$csize)
g_wheat <- induced_subgraph(g, which(comps$membership == biggest))

cat(sprintf("Largest component: %d nodes, %d edges\n", vcount(g_wheat), ecount(g_wheat)))

# ── Node attributes ────────────────────────────────────────────────────────────

shorten <- function(x) gsub("f\\. sp\\. ", "f.sp.", x)
V(g_wheat)$label <- shorten(V(g_wheat)$name)
V(g_wheat)$order <- nodes$order[match(V(g_wheat)$name, nodes$name)]

# ── Layout ─────────────────────────────────────────────────────────────────────

set.seed(42)
layout <- create_layout(g_wheat, layout = "stress")
layout <- stretch_center(layout, power = 0.55)

# ── Hull data: one row per (node, host) pair ───────────────────────────────────
# Collapse subspecies/variety to binomial so hulls group by species.

node_hosts_raw <- read.delim(NODE_HOSTS_TSV, stringsAsFactors = FALSE)
# Collapse subspecies/variety to binomial
node_hosts_raw$host <- sub("^(\\S+\\s+\\S+)\\s+.*", "\\1", node_hosts_raw$host)
# Remap Aegilops tauschii (D-genome donor) and all Triticum species to "Triticum sp."
# Aegilops STAT calls are indistinguishable from wheat D-subgenome reads.
# Collapsing Triticum species avoids spurious hull splits between T. aestivum / T. dicoccoides etc.
node_hosts_raw$host[grepl("^Aegilops", node_hosts_raw$host)] <- "Triticum sp."
node_hosts_raw$host[grepl("^Triticum", node_hosts_raw$host)] <- "Triticum sp."
node_hosts_raw$host[node_hosts_raw$host == "Triticinae"]     <- "Triticum sp."
# Re-aggregate after remapping
node_hosts_raw <- node_hosts_raw |>
  group_by(name, host) |>
  summarise(n = sum(n), .groups = "drop") |>
  as.data.frame()

cluster_names <- V(g_wheat)$name
layout_xy <- data.frame(name = layout$name, x = layout$x, y = layout$y,
                        stringsAsFactors = FALSE)

hull_data <- node_hosts_raw |>
  filter(name %in% cluster_names, n >= MIN_HOST_RUNS,
         host %in% HULL_WHITELIST) |>
  group_by(host) |>
  filter(n() >= MIN_HOST_NODES) |>
  ungroup() |>
  inner_join(layout_xy, by = "name")

present_hosts <- sort(unique(hull_data$host))
hull_cols     <- HOST_PALETTE[names(HOST_PALETTE) %in% present_hosts]
missing       <- setdiff(present_hosts, names(HOST_PALETTE))
if (length(missing) > 0) {
  hull_cols <- c(hull_cols, setNames(rep("#d8d8d8", length(missing)), missing))
}

cat(sprintf("Hull hosts (%d): %s\n", length(present_hosts),
            paste(present_hosts, collapse = ", ")))

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggraph(layout) +
  geom_mark_hull(
    aes(x, y, fill = host),
    data      = hull_data,
    concavity = 3,
    expand    = unit(6, "mm"),
    alpha     = 0.40
  ) +
  scale_fill_manual(values = hull_cols, name = "Host Plant",
                    guide = guide_legend(order = 4,
                      override.aes = list(alpha = 0.6))) +
  new_scale_fill() +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour      = "#555555",
    show.legend = c(width = TRUE, alpha = FALSE)
  ) +
  scale_edge_width(
    range  = c(0.3, 5),
    name   = "Co-occurrence Frequency",
    breaks = c(1, 5, 15, 30),
    guide  = guide_legend(order = 2)
  ) +
  scale_edge_alpha(range = c(0.20, 0.90), guide = "none") +
  geom_node_point(
    aes(size = total, fill = order),
    shape = 21, colour = "white", stroke = 0.4, alpha = 0.95
  ) +
  geom_node_label(
    aes(label = label, colour = order),
    size          = 4.5,
    label.padding = unit(0.20, "lines"),
    fill          = "white",
    alpha         = 0.88,
    repel         = TRUE,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +
  scale_fill_manual(values = ORDER_COLS, name = "Pathogen Order",
                    guide = guide_legend(order = 1,
                      override.aes = list(size = 4, shape = 21,
                                          stroke = 0.4))) +
  scale_colour_manual(values = ORDER_COLS, guide = "none") +
  scale_size_continuous(
    name   = "Pathogen Detections",
    range  = c(3, 14),
    breaks = c(5, 15, 30, 60, 80),
    guide  = guide_legend(order = 3,
               override.aes = list(shape = 21, fill = "grey55",
                                   colour = "white", stroke = 0.4, alpha = 0.95))
  ) +
  labs(
    title    = "Co-infection Network of Cereal Associated pathogens in Field RNA-seq data",
    subtitle = paste0(
      "Each edge connects pathogens co-detected in the same field biosample",
      " (NCBI SRA); edge width proportional to co-occurrence frequency;\n",
      "node size = detection count; node colour = fungal order;",
      " overlapping hulls = host plants spanning >= 3 pathogens in cluster (>= 2 runs each)"
    )
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 20, face = "bold"),
    plot.subtitle   = element_text(size = 13, colour = "grey40"),
    legend.text     = element_text(size = 13),
    legend.title    = element_text(size = 15, face = "bold"),
    legend.key.size = unit(1.1, "lines")
  )

ggsave(OUT_PDF, p, width = 16, height = 12)
ggsave(OUT_PNG, p, width = 16, height = 12, dpi = 200)
cat(sprintf("Written: %s\n", OUT_PNG))
