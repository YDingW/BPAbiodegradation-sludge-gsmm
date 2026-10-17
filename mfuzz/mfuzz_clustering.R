# Mfuzz soft clustering of the temporal profiles of the MAGs.
#
# Input (data/):
#   mag_abundance.csv   MAG x sampling point: metaWRAP abundance of the 90 MAGs (genome copies per million reads)
#   mag_expression.csv  MAG x sampling point: total transcript abundance of the 89 MAGs with transcripts (absolute quantification)
# Output (results/, <prefix> = abundance or expression):
#   <prefix>_clusters.csv      cluster of every MAG, its membership in every cluster, maximal membership, core (membership >= 0.5)
#   <prefix>_profiles.csv      standardised profile of every MAG with its cluster and membership (long table; input of trend_plot.R)
#   <prefix>_centroids.csv     cluster centroids (standardised units)
#   <prefix>_cluster_scan.csv  minimal distance between centroids and mean maximal membership for c = 2..10 (choice of c)
#   clusters_abundance_vs_expression.csv  both cluster assignments of every MAG side by side
# Method: log10(x + 1) -> standardisation of every MAG profile to mean 0 and sd 1 (Mfuzz::standardise)
#         -> fuzzy c-means (Mfuzz::mfuzz) with c = 3 clusters and the fuzzifier m of Mfuzz::mestimate;
#         the random start of the algorithm is fixed with set.seed(2026).
# Environment: R 4.5.2 with Biobase 2.70.0 and Mfuzz 2.70.0 (Bioconductor) and e1071 1.7-17.
# Run from this folder:  Rscript mfuzz_clustering.R

library(Biobase)
library(Mfuzz)

cluster_profiles <- function(file, prefix, c = 3, seed = 2026, core_cutoff = 0.5) {
  x <- as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
  eset <- standardise(ExpressionSet(assayData = log10(x + 1)))
  m <- mestimate(eset)
  # how c was chosen: distance between the closest centroids and mean maximal membership for c = 2..10
  scan <- t(sapply(2:10, function(k) {
    set.seed(seed)
    cl <- mfuzz(eset, c = k, m = m)
    d <- as.matrix(dist(cl$centers)); diag(d) <- NA
    c(c = k, min_centroid_distance = min(d, na.rm = TRUE), mean_max_membership = mean(apply(cl$membership, 1, max)))
  }))
  write.csv(scan, file.path("results", paste0(prefix, "_cluster_scan.csv")), row.names = FALSE)
  # the clustering
  set.seed(seed)
  cl <- mfuzz(eset, c = c, m = m)
  membership <- cl$membership
  colnames(membership) <- paste0("membership_", seq_len(c))
  clusters <- data.frame(MAG = rownames(x), cluster = cl$cluster, membership,
                         membership = apply(membership, 1, max), row.names = NULL, check.names = FALSE)
  clusters$core <- clusters$membership >= core_cutoff
  write.csv(clusters, file.path("results", paste0(prefix, "_clusters.csv")), row.names = FALSE)
  z <- exprs(eset)
  profiles <- data.frame(MAG = rep(rownames(z), ncol(z)), cluster = rep(cl$cluster, ncol(z)),
                         membership = rep(clusters$membership, ncol(z)), core = rep(clusters$core, ncol(z)),
                         time = rep(colnames(z), each = nrow(z)), zscore = as.vector(z))
  write.csv(profiles, file.path("results", paste0(prefix, "_profiles.csv")), row.names = FALSE)
  write.csv(data.frame(cluster = seq_len(c), cl$centers, check.names = FALSE), file.path("results", paste0(prefix, "_centroids.csv")),
            row.names = FALSE)
  cat(sprintf("%s: %d MAGs, %d sampling points, m = %.3f, c = %d; cluster sizes %s; core members %s\n", prefix, nrow(x), ncol(x), m, c,
              paste(table(cl$cluster), collapse = " / "), paste(tapply(clusters$core, clusters$cluster, sum), collapse = " / ")))
  clusters[, c("MAG", "cluster", "membership")]
}

a <- cluster_profiles("data/mag_abundance.csv", "abundance")
e <- cluster_profiles("data/mag_expression.csv", "expression")
both <- merge(a, e, by = "MAG", all = TRUE, suffixes = c("_abundance", "_expression"))
names(both) <- c("MAG", "abundance_cluster", "abundance_membership", "expression_cluster", "expression_membership")
write.csv(both, "results/clusters_abundance_vs_expression.csv", row.names = FALSE)
print(table(abundance = both$abundance_cluster, expression = both$expression_cluster))
