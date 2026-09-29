library(ggplot2)
d  <- read.delim("kraken/assign/figures/complexity_byhost.tsv", sep="\t", check.names=FALSE)
sl <- read.delim("kraken/assign/figures/model2d.tsv", sep="\t", check.names=FALSE)

# four host-flip cases: same organism real on one host, false on another
# clean textbook cases: two real (tight rising), one host-flip, one clean false
pick <- list(
  c("Zymoseptoria tritici","Triticum aestivum","Triticum aestivum"),
  c("Puccinia triticina","Triticum aestivum","Triticum aestivum"),
  c("Phakopsora pachyrhizi","Glycine max","Triticum aestivum"),
  c("Melampsora laricis-populina","Triticum aestivum","Triticum aestivum"))
keep <- do.call(rbind, lapply(pick, function(p)
  d[d$species==p[1] & d$host %in% p[2:3], ]))
keep <- keep[!duplicated(keep[c("species","host","reads","distinct")]), ]
keep$species <- factor(keep$species, levels=sapply(pick, `[`, 1))

PAL <- c("Triticum aestivum"="#eb6834","Glycine max"="#2a78d6","Zea mays"="#1baf7a",
         "Brassica napus"="#4a3aa7")
# per-host slope labels, stacked top-left
sl2 <- sl[paste(sl$species,sl$host) %in% unique(unlist(lapply(pick,function(p) paste(p[1],p[2:3])))), ]
sl2$species <- factor(sl2$species, levels=sapply(pick,`[`,1))
sl2 <- sl2[order(sl2$species, sl2$host), ]
sl2$row <- ave(seq_len(nrow(sl2)), sl2$species, FUN=seq_along)
sl2$lab <- sprintf("%s: slope %.2f, %d k-mers", sub(" .*","",sl2$host), sl2$slope, sl2$med_distinct)

p <- ggplot(keep, aes(reads, distinct, colour=host)) +
  geom_abline(slope=1, intercept=log10(0.03), linetype="dashed", linewidth=0.7, colour="grey55") +
  geom_point(size=1.1, alpha=0.45) +
  geom_text(data=sl2, aes(x=100, y=10^(7-0.6*(row-1)), label=lab, colour=host),
            hjust=0, size=3.5, fontface="bold", show.legend=FALSE) +
  facet_wrap(~species, ncol=2) +
  scale_x_log10(breaks=10^c(2,4,6), labels=c("100","10k","1M")) +
  scale_y_log10(breaks=10^c(1,4,7), labels=c("10","10k","10M")) +
  scale_colour_manual(values=PAL, name=NULL) +
  guides(colour=guide_legend(override.aes=list(size=3.5, alpha=1))) +
  labs(x="reads assigned to the organism", y="distinct genome fragments seen (k-mers)",
       title="Telling real infections from noise",
       subtitle="Rising cloud = real (genome coverage grows with more data). Flat cloud = false (reads stack on one fragment).\nPhakopsora is the same organism on two hosts: real on soybean, noise on wheat.") +
  theme_minimal(base_size=14) +
  theme(panel.background=element_rect(fill="white",colour=NA),
        plot.background=element_rect(fill="white",colour=NA),
        strip.text=element_text(size=12, face="italic"),
        legend.position="top", plot.title=element_text(size=15),
        plot.subtitle=element_text(size=11))
ggsave("true_detections/example_detections.png", p, width=12, height=9, dpi=200,
       bg="white", device=ragg::agg_png)
cat("ok\n")
