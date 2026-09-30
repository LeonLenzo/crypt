library(ggplot2)
d <- read.delim("kraken/filter/data/three.tsv", sep="\t", check.names=FALSE)
d$species <- factor(d$species, levels=unique(d$species))
d$method  <- factor(d$method, levels=c("A line","B cluster","C rarefaction"))
p <- ggplot(d, aes(reads, distinct, colour=cls)) +
  geom_point(size=0.8, alpha=0.4) +
  facet_grid(species ~ method, switch="y") +
  scale_x_log10(breaks=10^c(2,4,6), labels=c("100","10k","1M")) +
  scale_y_log10(breaks=10^c(1,3,5,7), labels=c("10","1k","100k","10M")) +
  scale_colour_manual(values=c(real="#2a78d6", `crash-out`="#eb6834"), name=NULL) +
  guides(colour=guide_legend(override.aes=list(size=4, alpha=1))) +
  labs(x="reads assigned to the organism", y="distinct genome fragments seen (k-mers)",
       title="Three ways to separate the real signal from the crash-out",
       subtitle="Rows = species (on wheat); columns = method. Blue = kept as real, orange = rejected.\nClustering (B) keeps 53% of the Phakopsora floor — no anchor for what 'real' is. A and C zero it correctly.") +
  theme_minimal(base_size=13) +
  theme(panel.background=element_rect(fill="white",colour=NA),
        plot.background=element_rect(fill="white",colour=NA),
        strip.text.y.left=element_text(size=9, face="italic", angle=90),
        strip.text.x=element_text(size=11),
        legend.position="top",
        plot.title=element_text(size=15), plot.subtitle=element_text(size=10.5))
ggsave("kraken/filter/figures/method_comparison.png", p, width=12, height=13, dpi=170,
       bg="white", device=ragg::agg_png)
cat("ok\n")
