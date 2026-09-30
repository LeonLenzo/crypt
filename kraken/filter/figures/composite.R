# Composite grid: one mini reads-vs-distinct chart per (pathogen x host) cell,
# rows and columns ordered by taxonomy, with the taxonomic trees on the margins.
# Pathogen tree (left) + host tree (top); each cell coloured by real-score.
suppressMessages({library(ggplot2); library(ape); library(ggtree); library(patchwork)})

d  <- read.delim("kraken/filter/data/_classified.tsv", sep="\t", check.names=FALSE)
pl <- read.delim("kraken/filter/data/patho_lineage.tsv", sep="\t", check.names=FALSE)
hl <- read.delim("kraken/filter/data/host_lineage.tsv",  sep="\t", check.names=FALSE)
for (c in c("score","reads","distinct")) d[[c]] <- suppressWarnings(as.numeric(d[[c]]))

cl <- d[d$cls!="unclassified" & is.finite(d$score), ]
keep_sp <- intersect(pl$label, unique(cl$species))
keep_ho <- intersect(hl$label, unique(cl$host))
pl <- pl[pl$label %in% keep_sp, ]; hl <- hl[hl$label %in% keep_ho, ]
d  <- d[d$species %in% keep_sp & d$host %in% keep_ho & d$reads>0 & d$distinct>0, ]

mktree <- function(x) {
  for (r in names(x)) x[[r]] <- factor(x[[r]])
  collapse.singles(as.phylo(~superkingdom/kingdom/phylum/class/order/family/genus/label,
                            data=x, collapse=FALSE))
}
tp <- mktree(pl); th <- mktree(hl)

# both trees built the same way (tips at y=1..n, x=depth); the pathogen tree is rotated
# into the top margin at draw time. Columns = pathogens, rows = hosts.
gp <- ggtree(tp, branch.length="none")
gh <- ggtree(th, branch.length="none")
px <- gp$data[gp$data$isTip, ]; sx <- setNames(px$y, px$label)   # pathogen column centres
hy <- gh$data[gh$data$isTip, ]; sy <- setNames(hy$y, hy$label)   # host row centres

d$sx <- sx[d$species]; d$sy <- sy[d$host]
nc <- length(sx); nr <- length(sy)
XL <- c(0.5, nc+0.5); YL <- c(0.5, nr+0.5)

resc <- function(v) pmin(pmax(log10(v)/8, 0), 1)                 # log10 value /8 -> [0,1]
d$gx <- d$sx + (resc(d$reads)-0.5)*0.9
d$gy <- d$sy + (resc(d$distinct)-0.5)*0.9
cells <- expand.grid(sx=sx, sy=sy)

abbr <- function(v) sub("^(\\S)\\S+ ", "\\1. ", v)

main <- ggplot() +
  geom_tile(data=cells, aes(sx, sy), fill="grey97", colour="grey88",
            linewidth=0.15, width=0.98, height=0.98) +
  geom_point(data=d, aes(gx, gy, colour=score), size=0.28, alpha=0.55) +
  scale_colour_gradient2(midpoint=0.5, low="#eb6834", mid="grey85", high="#2a78d6",
                         limits=c(0,1), name="real score", na.value="grey80") +
  scale_x_continuous(breaks=sx, labels=abbr(names(sx)), limits=XL, expand=c(0,0)) +
  scale_y_continuous(breaks=sy, labels=names(sy), limits=YL, expand=c(0,0), position="right") +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=8) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6, face="italic"),
        axis.text.y=element_text(size=6, face="italic"),
        legend.position="right",
        plot.margin=margin(2,2,2,2))

left <- gh + coord_cartesian(ylim=YL, expand=FALSE) + theme(plot.margin=margin(2,0,2,2))
# pathogen tree rotated into the top margin: tip index -> x, depth pointing down to the cells
top <- gp + coord_flip(ylim=XL, expand=FALSE) + scale_y_reverse() +
  theme(plot.margin=margin(2,2,0,2))

# explicit areas: plots map to these in provision order (top, guide, left, main)
layout <- c(area(t=1, l=2, b=1, r=2),   # host tree, top-centre
            area(t=1, l=3, b=2, r=3),   # legend, right column
            area(t=2, l=1, b=2, r=1),   # pathogen tree, left column
            area(t=2, l=2, b=2, r=2))   # matrix, centre
p <- top + guide_area() + left + main +
  plot_layout(design=layout, widths=c(1.4, 9, 1.0), heights=c(1.3, 9), guides="collect") +
  plot_annotation(
    title="Detection clouds by taxonomy: host (rows) x pathogen (columns)",
    subtitle="Each cell: log10 reads (x) vs log10 unique k-mers (y), axes 0-8. Blue rises (real), orange flat (crash-out).",
    theme=theme(plot.title=element_text(size=13), plot.subtitle=element_text(size=9)))

W <- 0.42*nc+4; H <- 0.36*nr+3
ggsave("kraken/filter/figures/composite_grid.png", p,
       width=W, height=H, dpi=150, limitsize=FALSE, bg="white", device=ragg::agg_png)
ggsave("kraken/filter/figures/composite_grid.pdf", p,          # vector, for zooming into cells
       width=W, height=H, limitsize=FALSE, bg="white", device=cairo_pdf)
cat("wrote composite_grid.{png,pdf}:", nr, "hosts x", nc, "pathogens\n")
