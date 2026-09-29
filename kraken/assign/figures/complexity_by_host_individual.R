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
  # one panel per host; strip carries the verdict so each is judged on its own
  ss$strip <- sprintf("%s\nR² %.2f · slope %.2f · %s",
                      ss$host, ss$r2, ss$slope,
                      ifelse(ss$verdict=="real (tight)", "✓ real", "✗ artefact"))
  ds <- merge(ds, ss[c("host","strip","verdict")], by = "host")
  ds$strip <- factor(ds$strip, levels = ss$strip[order(-ss$r2)])
  ncol <- min(3, length(hosts))

  p <- ggplot(ds, aes(reads, distinct, colour = verdict)) +
    geom_point(size = 1.3, alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, formula = y ~ x) +
    facet_wrap(~strip, ncol = ncol) +
    scale_x_log10(breaks = 10^c(2,4,6), labels = c("100","10k","1M")) +
    scale_y_log10(breaks = 10^c(1,3,5,7), labels = c("10","1k","100k","10M")) +
    scale_colour_manual(values = c(`real (tight)`="#2a78d6", `artefact (loose)`="#eb6834"),
                        guide = "none") +
    labs(x = "reads assigned to the organism", y = "distinct genome fragments seen (k-mers)",
         title = bquote(italic(.(sp))),
         subtitle = "One panel per host. Tight rising fit (blue) = real infection; flat or scattered (red) = artefact.") +
    theme_minimal(base_size = 14) +
    theme(panel.background = element_rect(fill="white", colour=NA),
          plot.background = element_rect(fill="white", colour=NA),
          strip.text = element_text(size = 10, lineheight = 1.05),
          plot.title = element_text(size = 17),
          plot.subtitle = element_text(size = 11))

  nrow <- ceiling(length(hosts)/ncol)
  fn <- gsub("[^A-Za-z0-9]+", "_", sp)
  ggsave(sprintf("complexity_by_host/%s.png", fn), p,
         width = 3.6*ncol + 1, height = 3.3*nrow + 1.2,
         dpi = 150, bg = "white", device = ragg::agg_png, limitsize = FALSE)
}
cat("wrote", length(unique(st$species)), "images to complexity_by_host/\n")
