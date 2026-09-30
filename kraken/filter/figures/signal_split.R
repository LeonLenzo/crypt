library(ggplot2)
d<-read.delim("kraken/filter/data/splitfig.tsv",sep="\t",check.names=FALSE)
d$panel<-factor(d$panel,levels=unique(d$panel))
lines<-unique(d[c("panel","a","b")])
p<-ggplot(d,aes(reads,distinct))+
  geom_point(aes(colour=cls),size=1.3,alpha=0.45)+
  geom_abline(data=lines,aes(slope=b,intercept=a),linetype="dashed",linewidth=0.9,colour="grey25")+
  facet_wrap(~panel,ncol=2)+
  scale_x_log10(breaks=10^c(2,4,6),labels=c("100","10k","1M"))+
  scale_y_log10(breaks=10^c(1,3,5,7),labels=c("10","1k","100k","10M"))+
  scale_colour_manual(values=c(real="#2a78d6",`crash-out`="#eb6834"),name=NULL)+
  guides(colour=guide_legend(override.aes=list(size=4,alpha=1)))+
  labs(x="reads assigned to the organism",y="distinct genome fragments seen (k-mers)",
       title="Separating the real signal from the crash-out, per detection",
       subtitle="Dashed = fitted real accumulation line. Blue = on the line (real); orange = fallen to the floor (crash-out).\nKeeps the real subset of a sporadic pathogen; a pure floor (Phakopsora/wheat) keeps nothing.")+
  theme_minimal(base_size=14)+
  theme(panel.background=element_rect(fill="white",colour=NA),plot.background=element_rect(fill="white",colour=NA),
        strip.text=element_text(size=12,face="italic"),legend.position="top",
        plot.title=element_text(size=15),plot.subtitle=element_text(size=10.5))
ggsave("kraken/filter/figures/signal_vs_crashout.png",p,width=12,height=9,dpi=200,bg="white",device=ragg::agg_png)
cat("ok\n")
