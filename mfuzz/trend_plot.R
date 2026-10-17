# Trend plot of the Mfuzz clusters: one panel per cluster, one line per MAG coloured by its membership,
# the black line is the cluster centroid of Mfuzz (membership-weighted mean profile).
# Input:  results/<prefix>_profiles.csv and results/<prefix>_centroids.csv, <prefix> = abundance or expression
#         (from mfuzz_clustering.R)
# Output: results/trend_plot.png, 600 dpi (abundance clusters AC1-3 on the left, expression clusters EC1-3 on the right)
# Environment: R 4.5.2 with ggplot2 4.0.3, ggh4x 0.3.1 and patchwork 1.3.2.
# Run from this folder:  Rscript trend_plot.R

library(ggplot2)
library(patchwork)

strip_colours <- list(abundance = c("#A6CEE3", "#B2DF8A", "#CAB2D6"), expression = c("#FDBF6F", "#FB9A99", "#D9D9D9"))
time_label <- function(x) gsub("([0-9])([A-Za-z])", "\\1 \\2", x)      # "24h" -> "24 h", "0D" -> "0 D"

trend_panel <- function(prefix, label, title) {
  d <- read.csv(file.path("results", paste0(prefix, "_profiles.csv")))
  times <- unique(d$time)
  d$time <- factor(time_label(d$time), levels = time_label(times))
  d$cluster <- factor(paste0(label, d$cluster))
  d <- d[order(d$membership), ]                       # low membership drawn first, high membership on top
  cen <- read.csv(file.path("results", paste0(prefix, "_centroids.csv")), check.names = FALSE)
  cen <- data.frame(cluster = factor(paste0(label, rep(cen$cluster, length(times)))),
                    time = factor(time_label(rep(times, each = nrow(cen))), levels = levels(d$time)),
                    zscore = unlist(cen[times]))
  strips <- ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = strip_colours[[prefix]]),
                                text_y = ggh4x::elem_list_text(face = "bold", size = 15))
  ggplot(d, aes(time, zscore, group = MAG)) +
    geom_line(aes(colour = membership, alpha = core), linewidth = 0.4) +
    scale_colour_viridis_c(option = "C", limits = c(0, 1), name = "Membership") +
    scale_alpha_manual(values = c("TRUE" = 0.95, "FALSE" = 0.25), guide = "none") +
    geom_line(data = cen, aes(time, zscore, group = cluster), colour = "black", linewidth = 1.4) +
    geom_point(data = cen, aes(time, zscore, group = cluster), colour = "black", size = 2.2) +
    ggh4x::facet_wrap2(~ cluster, ncol = 1, strip.position = "left", strip = strips) +
    scale_y_continuous(position = "right") +
    labs(title = title, x = NULL, y = "Standardized abundance / expression") +
    theme_bw(base_size = 16, base_family = "Arial") +
    theme(text = element_text(face = "bold", colour = "black"), axis.text = element_text(colour = "black", size = 13),
          axis.title = element_text(size = 15), plot.title = element_text(size = 18, hjust = 0),
          legend.text = element_text(size = 13), legend.title = element_text(size = 14),
          panel.border = element_rect(colour = "grey50", linewidth = 0.5), panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_blank())
}

p <- trend_panel("abundance", "AC", "Abundance") | trend_panel("expression", "EC", "Expression")
ggsave("results/trend_plot.png", p, width = 10, height = 6.4, units = "in", dpi = 600, bg = "white")
