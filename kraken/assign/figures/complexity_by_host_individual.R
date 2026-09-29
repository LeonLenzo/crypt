# complexity_by_host_individual.R — one image per species, reflecting the tightness model.
#
# reads vs distinct minimizers, points + fitted line per host. A host where the fit is
# tight (R2 >= 0.5, slope > 0) is a real infection (blue line); a loose or flat fit is a
# pile-up or cross-map (red line). Tightness, not level, so low-abundance real detections
# are kept. One PNG per species in complexity_by_host/.

library(ggplot2)
d  <- read.delim("complexity_by_host/_points.tsv", sep="\t", check.names=FALSE)
st <- read.delim("complexity_by_host/_stats.tsv",  sep="\t", check.names=FALSE)
R_CUT <- 0.5
st$verdict <- ifelse(st$r2 >= R_CUT & st$slope > 0, "real (tight)", "artefact (loose)")

HOSTPAL <- c("#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#4a3aa7","#008300",
             "#00a0b0","#a05a2c","#777777")

for (sp in sort(unique(st$species))) {          # only species with a fittable host
  ds <- d[d$species == sp, ]
  ss <- st[st$species == sp, ]
  # keep hosts that have a fitted stat; fold the rest out
  ds <- ds[ds$host %in% ss$host, ]
  hosts <- ss$host[order(-ss$r2)]
  pal <- setNames(HOSTPAL[seq_along(hosts)], hosts)
  ds$host <- factor(ds$host, levels = hosts)
  ss$host <- factor(ss$host, levels = hosts)
  ss <- ss[order(ss$host), ]
  ss$lab <- sprintf("%s  —  R² %.2f, slope %.2f  %s",
                    ss$host, ss$r2, ss$slope,
                    ifelse(ss$verdict=="real (tight)", "✓ real", "✗ artefact"))
  ss$y <- 10^(max(log10(ds$distinct)) * (1 - 0.07*(seq_len(nrow(ss))-1)))

  p <- ggplot(ds, aes(reads, distinct, colour = host)) +
    geom_point(size = 1.4, alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, formula = y ~ x) +
    geom_text(data = ss, aes(x = min(ds$reads), y = y, label = lab, colour = host),
              hjust = 0, size = 4, fontface = "bold", show.legend = FALSE) +
    scale_x_log10(breaks = 10^(2:7), labels = c("100","1k","10k","100k","1M","10M")) +
    scale_y_log10(breaks = 10^(1:7), labels = c("10","100","1k","10k","100k","1M","10M")) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = "reads assigned to the organism", y = "distinct genome fragments seen (k-mers)",
         title = bquote(italic(.(sp))),
         subtitle = "Tight rising fit = real infection. Flat or scattered fit = artefact. Coverage tracks depth, not abundance.") +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill="white", colour=NA),
          plot.background = element_rect(fill="white", colour=NA),
          plot.title = element_text(size = 17),
          plot.subtitle = element_text(size = 11))

  fn <- gsub("[^A-Za-z0-9]+", "_", sp)
  ggsave(sprintf("complexity_by_host/%s.png", fn), p, width = 10, height = 7,
         dpi = 150, bg = "white", device = ragg::agg_png)
}
cat("wrote", length(unique(st$species)), "images to complexity_by_host/\n")
