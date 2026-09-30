# Co-occurrence networks, one per host, panelled together.
# Node = pathogen (size = biosample prevalence, colour = taxonomic class);
# edge = co-occurrence in a biosample (width = number of shared biosamples).
suppressMessages({library(igraph); library(ggraph); library(ggplot2); library(patchwork)})
set.seed(1)
PADJ <- 0.05                                           # keep only associations beyond chance

nd <- read.delim("kraken/filter/data/network_nodes.tsv", sep="\t", check.names=FALSE)
ed <- read.delim("kraken/filter/data/network_edges.tsv", sep="\t", check.names=FALSE)
sm <- read.delim("kraken/filter/data/network_summary.tsv", sep="\t", check.names=FALSE)
ed <- ed[ed$p_adj < PADJ & ed$log2_oe > 0, ]           # significant POSITIVE associations
ed$weight <- ed$log2_oe                                # edge width = fold-enrichment over chance

# class palette pinned to class names, so a class is always the same colour
pal <- c(Agaricomycetes="#4477AA", Dothideomycetes="#EE6677", Eurotiomycetes="#228833",
         Leotiomycetes="#CCBB44", Peronosporomycetes="#66CCEE", Pucciniomycetes="#AA3377",
         Sordariomycetes="#BBBBBB", Ustilaginomycetes="#EE8866", other="#000000")
cls <- names(pal)[names(pal) %in% nd$class]
abbr <- function(v) sub("^(\\S)\\S+ ", "\\1. ", v)

hosts <- sm$host[order(-sm$bs_coinfected)]
plots <- list()
for (h in hosts) {
  e <- ed[ed$host==h, c("source","target","weight")]
  if (nrow(e) < 3) next                                # skip hosts with almost no structure
  n <- nd[nd$host==h, ]
  keep <- union(e$source, e$target)                   # drop isolated nodes
  n <- n[n$species %in% keep, ]
  n$class <- factor(n$class, levels=cls)               # identical scale across panels -> one legend
  g <- graph_from_data_frame(e, vertices=n[,c("species","n_biosamples","class")], directed=FALSE)
  denom <- sm$bs_total[sm$host==h]; nco <- sm$bs_coinfected[sm$host==h]
  p <- ggraph(g, layout="stress") +
    geom_edge_link(aes(width=weight), colour="grey70", alpha=0.5, show.legend=FALSE) +
    geom_node_point(aes(size=n_biosamples, fill=class), shape=21, colour="grey30", stroke=0.3) +
    geom_node_text(aes(label=abbr(name)), size=2.1, repel=TRUE, max.overlaps=Inf,
                   segment.size=0.2, colour="grey20") +
    scale_edge_width(range=c(0.4, 3), guide="none") +
    scale_size_area(max_size=10, guide="none") +
    scale_fill_manual(values=pal, limits=cls, drop=FALSE, name="pathogen class") +
    guides(fill=guide_legend(override.aes=list(size=4))) +
    labs(title=sprintf("%s", h),
         subtitle=sprintf("%d biosamples, %d co-infected · edge = co-occurrence beyond chance (Fisher BH p < %.2f)", denom, nco, PADJ)) +
    theme_void(base_size=12) +
    theme(plot.title=element_text(size=15, face="italic"),
          plot.subtitle=element_text(size=9, colour="grey40"),
          legend.position="right", plot.margin=margin(8,8,8,8))
  plots[[h]] <- p
}

# one image per host
outdir <- "kraken/filter/figures/networks"
dir.create(outdir, showWarnings=FALSE)
invisible(file.remove(list.files(outdir, pattern="\\.png$", full.names=TRUE)))  # drop stale hosts
for (h in names(plots)) {
  f <- file.path(outdir, gsub("[^A-Za-z0-9]+", "_", h))
  ggsave(paste0(f, ".png"), plots[[h]], width=9, height=6.5, dpi=150,
         bg="white", device=ragg::agg_png)
}
cat("wrote", length(plots), "host networks to", outdir, "\n")
