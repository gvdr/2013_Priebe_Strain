theme_white_hc <- theme(panel.background = element_rect(fill = "white"),
                        panel.grid.minor = element_line(size = .3, color = "lightgrey"),
                        panel.grid.major = element_line(size = .5, color = "lightgrey"),
                        axis.line = element_line(size=.7, color = "grey"),
                        text = element_text(size=14),
                        strip.background = element_rect(fill= NA, color = NA),
                        strip.text = element_text(size=14, family = "Helvetica"))

#Hat tip to boB Rudis
#viridis_scales.R
#https://gist.github.com/hrbrmstr
viridis_pal <- function(alpha=1) {
  function(n) {
    viridis(n, alpha)
  }
}
scale_color_viridis <- function(..., alpha=1, discrete=TRUE) {
  if (discrete) {
    discrete_scale("colour", "viridis", viridis_pal(alpha), ...)
  } else {
    scale_color_gradientn(colours = viridis(256, alpha), ...)
  }
}
scale_fill_viridis <- function (..., alpha=1, discrete=TRUE) {
  if (discrete) {
    discrete_scale("fill", "viridis", viridis_pal(alpha), ...)
  } else {
    scale_fill_gradientn(colours = viridis(256, alpha), ...)
  }
}


ExportPlot <- function(gplot, filename, width=2, height=1.5) {
  # Export plot in PDF and EPS.
  # Notice that A4: width=11.69, height=8.27
  pdf(paste(filename,'.pdf', sep = ""),width = width, height = height,useDingbats = FALSE)
  plot(gplot)
  dev.off()
  ggsave(paste(filename, '.eps', sep=""), gplot)
  ggsave(paste(filename, '.svg', sep=""), gplot)
  ggsave(paste(filename, '.png', sep=""), gplot, width = width, height = height)
}

#Hat tip to Hadley
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  arrangeGrob(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}