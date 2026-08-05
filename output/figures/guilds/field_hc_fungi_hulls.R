#!/usr/bin/env Rscript
# figure/guilds/field_hc_fungi_hulls.R
# Fungi-only co-infection network with host hulls, high-confidence field co-infections
# (diff-genus secondaries only, same_genus_secondary=False). Run from crypt/.

library(igraph)
library(ggraph)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(dplyr)
library(graphlayouts)

NODES_TSV      <- "figure/guilds/field_hc_nodes.tsv"
EDGES_TSV      <- "figure/guilds/field_hc_edges.tsv"
NODE_HOSTS_TSV <- "figure/guilds/field_hc_node_hosts.tsv"
OUT_PDF        <- "figure/guilds/field_hc_fungi_hulls.pdf"
OUT_PNG        <- "figure/guilds/field_hc_fungi_hulls.png"

MIN_HOST_RUNS  <- 1   # min runs for a (pathogen, host) pair to count
MIN_HOST_NODES <- 3   # min distinct cluster nodes a host must span

HULL_WHITELIST <- c(
  "Triticum sp.", "Zea mays", "Hordeum vulgare", "Sorghum bicolor",
  "Vitis vinifera", "Theobroma cacao", "Solanum tuberosum",
  "Oryza sativa", "Cicer arietinum", "Mangifera indica"
)

# ── Taxonomy ──────────────────────────────────────────────────────────────────

# Warm palette for fungal class nodes
CLASS_COLS <- c(
  Dothideomycetes   = "#e67e22",
  Sordariomycetes   = "#c0392b",
  Leotiomycetes     = "#f1c40f",
  Eurotiomycetes    = "#a04000",
  Pucciniomycetes   = "#e91e63",
  Agaricomycetes    = "#784212",
  Ustilaginomycetes = "#ff7043",
  Unknown           = "#95a5a6"
)

# Cool palette for host plant hull borders
HOST_BORDER_COLS <- c(
  "Triticum sp."      = "#1565c0",
  "Zea mays"          = "#2e7d32",
  "Hordeum vulgare"   = "#00838f",
  "Sorghum bicolor"   = "#4527a0",
  "Vitis vinifera"    = "#6a1b9a",
  "Theobroma cacao"   = "#00695c",
  "Solanum tuberosum" = "#0277bd",
  "Oryza sativa"      = "#01579b",
  "Cicer arietinum"   = "#283593",
  "Mangifera indica"  = "#00796b"
)

CLASS_MAP <- c(
  # Pucciniomycetes (rusts)
  "Puccinia graminis"                    = "Pucciniomycetes",
  "Puccinia graminis f. sp. tritici"    = "Pucciniomycetes",
  "Puccinia striiformis"                = "Pucciniomycetes",
  "Puccinia striiformis f. sp. tritici" = "Pucciniomycetes",
  "Puccinia triticina"                  = "Pucciniomycetes",
  # Dothideomycetes — Pleosporales
  "Alternaria alternata"                = "Dothideomycetes",
  "Alternaria solani"                   = "Dothideomycetes",
  "Bipolaris maydis"                    = "Dothideomycetes",
  "Bipolaris sorokiniana"               = "Dothideomycetes",
  "Bipolaris zeicola"                   = "Dothideomycetes",
  "Cercospora beticola"                 = "Dothideomycetes",
  "Cercospora kikuchii"                 = "Dothideomycetes",
  "Cercospora zeina"                    = "Dothideomycetes",
  "Exserohilum turcica"                 = "Dothideomycetes",
  "Exserohilum turcicum"                = "Dothideomycetes",
  "Leptosphaeria maculans"              = "Dothideomycetes",
  "Parastagonospora nodorum"            = "Dothideomycetes",
  "Pyrenophora tritici-repentis"        = "Dothideomycetes",
  "Setosphaeria turcica"                = "Dothideomycetes",
  "Stemphylium lycopersici"             = "Dothideomycetes",
  # Dothideomycetes — Capnodiales
  "Zymoseptoria tritici"                = "Dothideomycetes",
  # Dothideomycetes — Botryosphaeriales
  "Lasiodiplodia theobromae"            = "Dothideomycetes",
  "Neofusicoccum parvum"                = "Dothideomycetes",
  # Sordariomycetes — Hypocreales
  "Fusarium asiaticum"                  = "Sordariomycetes",
  "Fusarium fujikuroi"                  = "Sordariomycetes",
  "Fusarium graminearum"                = "Sordariomycetes",
  "Fusarium odoratissimum"              = "Sordariomycetes",
  "Fusarium oxysporum"                  = "Sordariomycetes",
  "Fusarium oxysporum f.sp.lycopersici" = "Sordariomycetes",
  "Fusarium proliferatum"               = "Sordariomycetes",
  "Fusarium pseudograminearum"          = "Sordariomycetes",
  "Fusarium solani"                     = "Sordariomycetes",
  "Fusarium verticillioides"            = "Sordariomycetes",
  "Trichoderma virens"                  = "Sordariomycetes",
  # Sordariomycetes — Glomerellales
  "Colletotrichum fructicola"           = "Sordariomycetes",
  "Colletotrichum gloeosporioides"      = "Sordariomycetes",
  "Colletotrichum higginsianum"         = "Sordariomycetes",
  "Colletotrichum siamense"             = "Sordariomycetes",
  # Sordariomycetes — Magnaporthales
  "Pyricularia oryzae"                  = "Sordariomycetes",
  # Leotiomycetes — Helotiales
  "Botrytis cinerea"                    = "Leotiomycetes",
  "Sclerotinia sclerotiorum"            = "Leotiomycetes",
  # Eurotiomycetes — Eurotiales
  "Aspergillus flavus"                  = "Eurotiomycetes",
  "Aspergillus niger"                   = "Eurotiomycetes",
  "Penicillium digitatum"               = "Eurotiomycetes",
  "Penicillium expansum"                = "Eurotiomycetes",
  # Agaricomycetes — Agaricales
  "Moniliophthora perniciosa"           = "Agaricomycetes",
  # Ustilaginomycetes — Ustilaginales
  "Ustilago maydis"                     = "Ustilaginomycetes"
)

# ── Load ──────────────────────────────────────────────────────────────────────

nodes <- read.delim(NODES_TSV, stringsAsFactors = FALSE) %>%
  filter(kingdom == "Fungi")
fungi_names <- nodes$name
edges <- read.delim(EDGES_TSV, stringsAsFactors = FALSE) %>%
  filter(node1 %in% fungi_names, node2 %in% fungi_names)
nodes <- nodes %>% filter(name %in% unique(c(edges$node1, edges$node2)))

nodes$class <- CLASS_MAP[nodes$name]
nodes$class[is.na(nodes$class)] <- "Unknown"

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

comps   <- components(g)
biggest <- which.max(comps$csize)
g_main  <- induced_subgraph(g, which(comps$membership == biggest))
cat(sprintf("Main component: %d nodes, %d edges\n", vcount(g_main), ecount(g_main)))

shorten <- function(x) gsub("f\\. sp\\. ", "f.sp.", x)
V(g_main)$label <- shorten(V(g_main)$name)
V(g_main)$class <- nodes$class[match(V(g_main)$name, nodes$name)]

# ── Hull node sets (precomputed for seed optimisation) ────────────────────────

node_hosts_raw <- read.delim(NODE_HOSTS_TSV, stringsAsFactors = FALSE)
node_hosts_raw$host <- sub("^(\\S+\\s+\\S+)\\s+.*", "\\1", node_hosts_raw$host)
node_hosts_raw$host[grepl("^Aegilops|^Triticum|^Triticinae", node_hosts_raw$host)] <- "Triticum sp."
node_hosts_raw <- node_hosts_raw %>%
  group_by(name, host) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  as.data.frame()

cluster_names <- V(g_main)$name

hull_node_sets <- node_hosts_raw %>%
  filter(name %in% cluster_names, n >= MIN_HOST_RUNS, host %in% HULL_WHITELIST) %>%
  group_by(host) %>%
  filter(n() >= MIN_HOST_NODES) %>%
  ungroup() %>%
  select(name, host)

# ── Layout (seed optimised for hull compactness) ───────────────────────────────

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

best_seed  <- 42
best_score <- Inf
for (s in c(1, 7, 13, 21, 42, 77, 99, 123, 200, 314, 500, 999)) {
  set.seed(s)
  ly    <- create_layout(g_main, layout = "stress")
  ly    <- stretch_center(ly, power = 0.55)
  ly_xy <- data.frame(name = ly$name, x = ly$x, y = ly$y)
  score <- hull_node_sets %>%
    inner_join(ly_xy, by = "name") %>%
    group_by(host) %>%
    summarise(spread = sqrt(var(x) + var(y)), .groups = "drop") %>%
    pull(spread) %>%
    mean(na.rm = TRUE)
  if (score < best_score) { best_score <- score; best_seed <- s }
}
cat(sprintf("Best seed: %d (hull compactness: %.3f)\n", best_seed, best_score))
set.seed(best_seed)
layout <- create_layout(g_main, layout = "stress")
layout <- stretch_center(layout, power = 0.55)

# ── Hull data ─────────────────────────────────────────────────────────────────

layout_xy <- data.frame(name = layout$name, x = layout$x, y = layout$y,
                        stringsAsFactors = FALSE)

hull_data <- node_hosts_raw %>%
  filter(name %in% cluster_names, n >= MIN_HOST_RUNS, host %in% HULL_WHITELIST) %>%
  group_by(host) %>%
  filter(n() >= MIN_HOST_NODES) %>%
  ungroup() %>%
  inner_join(layout_xy, by = "name")

present_hosts <- sort(unique(hull_data$host))
hull_cols     <- HOST_BORDER_COLS[names(HOST_BORDER_COLS) %in% present_hosts]
missing       <- setdiff(present_hosts, names(HOST_BORDER_COLS))
if (length(missing) > 0)
  hull_cols <- c(hull_cols, setNames(rep("#999999", length(missing)), missing))

cat(sprintf("Hull hosts (%d): %s\n", length(present_hosts),
            paste(present_hosts, collapse = ", ")))

# ── Plot ──────────────────────────────────────────────────────────────────────

p <- ggraph(layout) +
  geom_mark_hull(
    aes(x, y, colour = host),
    data      = hull_data,
    fill      = NA,
    linewidth = 1.2,
    concavity = 3,
    expand    = unit(6, "mm")
  ) +
  scale_colour_manual(values = hull_cols, name = "Host Plant",
                      guide = guide_legend(order = 4,
                        override.aes = list(fill = NA, linewidth = 1.2))) +
  new_scale_colour() +
  geom_edge_link(
    aes(width = weight, alpha = weight),
    colour = "#555555",
    show.legend = c(width = TRUE, alpha = FALSE)
  ) +
  scale_edge_width(range = c(0.3, 5), name = "Co-occurrence Frequency",
                   breaks = c(1, 5, 15, 30),
                   guide = guide_legend(order = 2)) +
  scale_edge_alpha(range = c(0.20, 0.90), guide = "none") +
  geom_node_point(
    aes(size = total, fill = class),
    shape = 21, colour = "white", stroke = 0.4, alpha = 0.95
  ) +
  geom_node_label(
    aes(label = label, colour = class),
    size = 4.2, label.padding = unit(0.18, "lines"),
    fill = "white", alpha = 0.88,
    repel = TRUE, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = CLASS_COLS, name = "Fungal Class",
                    guide = guide_legend(order = 1,
                      override.aes = list(size = 4, shape = 21,
                                          stroke = 0.4))) +
  scale_colour_manual(values = CLASS_COLS, guide = "none") +
  scale_size_continuous(
    name = "Pathogen Detections", range = c(3, 14),
    breaks = c(5, 15, 30, 60, 80),
    guide = guide_legend(order = 3,
      override.aes = list(shape = 21, fill = "grey55",
                          colour = "white", stroke = 0.4, alpha = 0.95))
  ) +
  labs(
    title    = "Fungal co-infection guilds - field samples (high confidence)",
    subtitle = paste0(
      "Nodes = fungal pathogens; edges = co-occurrence in same field biosample (diff-genus secondaries only)\n",
      "Node colour = fungal class; hulls = host plants spanning >= 3 nodes"
    )
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 18, face = "bold"),
    plot.subtitle   = element_text(size = 11, colour = "grey40"),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 13, face = "bold"),
    legend.key.size = unit(1.0, "lines")
  )

ggsave(OUT_PDF, p, width = 16, height = 12)
ggsave(OUT_PNG, p, width = 16, height = 12, dpi = 200)
cat(sprintf("Written: %s\n", OUT_PNG))
