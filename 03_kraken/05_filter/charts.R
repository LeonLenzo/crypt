# One image per species: reads vs distinct k-mers, faceted by host, each detection
# labelled real (on the accumulation line) or crash-out (fallen to the floor).
library(ggplot2)
d <- read.delim("03_kraken/05_filter/data/_classified.tsv", sep="\t", check.names=FALSE)
pow10 <- function(x) round(log10(x))                       # axis labels: log10 exponent
LIM <- 10^c(0, 8)                                          # common limits across every chart

for (sp in sort(unique(d$species))) {
  ds <- d[d$species==sp, ]
  # only hosts we could actually classify (drop sparse 'too few' panels)
  classifiable <- names(which(tapply(ds$cls!="unclassified", ds$host, any)))
  ds <- ds[ds$host %in% classifiable, ]
  if (nrow(ds)==0) next
  # host panel labels with real fraction
  tab <- tapply(ds$cls, ds$host, function(v){
    n<-length(v); r<-sum(v=="real"); u<-sum(v=="unclassified")
    if (u==n) sprintf("%%s\n%d runs (too few to classify)", n)
    else sprintf("%%s\n%d runs · %d%%%% real", n, round(100*r/n))
  })
  ds$panel <- factor(ds$host, levels=names(sort(tapply(ds$cls=="real",ds$host,mean),decreasing=TRUE)))
  labs <- setNames(mapply(function(h) sprintf(tab[[h]], h), levels(ds$panel)), levels(ds$panel))
  ds$panel <- factor(labs[as.character(ds$host)], levels=labs)
  ncol <- min(3, length(unique(ds$host)))
  ds$score <- suppressWarnings(as.numeric(ds$score))
  p <- ggplot(ds, aes(reads, distinct, colour=score)) +
    geom_point(size=1.3, alpha=0.55) +
    facet_wrap(~panel, ncol=ncol) +
    scale_x_log10(breaks=10^c(0,4,8), labels=pow10, limits=LIM) +
    scale_y_log10(breaks=10^c(0,4,8), labels=pow10, limits=LIM) +
    scale_colour_gradient2(midpoint=0.5, low="#eb6834", mid="grey85", high="#2a78d6",
                           limits=c(0,1), name="real score", na.value="grey80") +
    labs(x="log10 reads assigned", y="log10 unique kmers",
         title=bquote(italic(.(sp))),
         subtitle="Weighted score (tightness+slope+coverage+on-line+unbiased). Blue = real (>=0.5), orange = crash-out.") +
    theme_minimal(base_size=13) +
    theme(panel.background=element_rect(fill="white",colour=NA),
          plot.background=element_rect(fill="white",colour=NA),
          strip.text=element_text(size=9, lineheight=1.05),
          legend.position="top", plot.title=element_text(size=16),
          plot.subtitle=element_text(size=10))
  nrow <- ceiling(length(unique(ds$host))/ncol)
  ggsave(sprintf("03_kraken/05_filter/charts/%s.png", gsub("[^A-Za-z0-9]+","_",sp)), p,
         width=3.7*ncol+1, height=3.2*nrow+1.1, dpi=150, bg="white",
         device=ragg::agg_png, limitsize=FALSE)
}
cat("wrote", length(unique(d$species)), "images\n")
