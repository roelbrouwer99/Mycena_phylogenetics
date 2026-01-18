# =========================
# README
# =========================
# Purpose: End-to-end tree-space report using RF distances + Ward.D2 clustering.
# Inputs:
#   - Two MrBayes posterior .t files (run1_path, run2_path)
#   - One ML tree for comparison (ml_tree_path)
# Outputs:
#   - 01_silhouette_scores_vs_K.pdf
#   - 02_MDS_colored_by_run.pdf
#   - 03_MDS_colored_by_cluster.pdf
#   - 04_dendrogram_by_cluster.pdf
#   - 05_RF_heatmap_ML_vs_clusterConsensus.pdf
#   - consensus_cluster_<i>_K<k>.tre and RF tables
# How to run:
#   Rscript tree_space.R
# Key settings to tweak:
#   burnin_frac, thin_every, K_min/K_max, consensus_p, input paths, outdir, base_family.
# Notes:
#   - Uses native R pipes (|>); requires R 4.1+.
# =========================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(cluster)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(dendextend)
})

# =========================
# 0) USER SETTINGS
# =========================

run1_path <- "C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/PRANK_Trees/BAYES/Perry_ITS/Perry_ITS_PRANK.best.fas.clipkit.cleaned.nex.run1.t"
run2_path <- "C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/PRANK_Trees/BAYES/Perry_ITS/Perry_ITS_PRANK.best.fas.clipkit.cleaned.nex.run2.t"

# Burn-in fraction (set to 0 to disable burn-in)
burnin_frac <- 0.25

# Combat speed constraints
thin_every <- 10

# K-search range for silhouette-based PAM clustering
K_min <- 2
K_max <- 10

# Consensus trees per cluster (at optimal K)
consensus_p <- 0.5

# ML comparison (ML vs cluster consensus RF)
ml_tree_path <- "C:/Users/rbrou/Desktop/Roel/Naturalis/Phylogenetics/3_Tree files/PRANK_Trees/IQ-TREE/contree/Perry_ITS_PRANK.best.fas.clipkit.contree"

# Output directory
outdir <- "tree_space_outputs_Perry_ITS_v3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- Shared plot style ---
plot_width <- 7
plot_height <- 6
plot_dpi <- 300

base_size <- 13
base_family <- "" # e.g., "Arial" (leave "" to use system default)

run_cols <- c("run1" = "firebrick3", "run2" = "dodgerblue3")

theme_report <- function() {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(p, filename, w = plot_width, h = plot_height) {
  ggsave(file.path(outdir, filename), p, width = w, height = h, dpi = plot_dpi)
}

# =========================
# 1) HELPERS
# =========================

normalize_labels <- function(x) {
  x <- trimws(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("-", "_", x)
  x <- gsub("_+", "_", x)
  toupper(x)
}

drop_fake_numeric_tips <- function(tr) {
  is_fake_tip <- grepl("^[0-9]+$", tr$tip.label)
  if (any(is_fake_tip)) {
    message("Removing ", sum(is_fake_tip), " fake numeric tips.")
    tr <- drop.tip(tr, tr$tip.label[is_fake_tip])
  }
  tr
}

make_cluster_palette <- function(K) {
  cols <- scales::hue_pal()(max(K, 2))
  names(cols) <- as.character(seq_len(max(K, 2)))
  cols
}

# Safe keep.tip that preserves order
keep_common_tips <- function(tr, common_tips) {
  keep.tip(tr, common_tips)
}

# =========================
# 2) LOAD TREES
# =========================

message("Loading MrBayes runs...")
run1 <- read.nexus(run1_path)
run2 <- read.nexus(run2_path)
message("Loaded trees: run1=", length(run1), ", run2=", length(run2))

burnin1 <- floor(burnin_frac * length(run1))
burnin2 <- floor(burnin_frac * length(run2))

post1 <- run1[(burnin1 + 1):length(run1)]
post2 <- run2[(burnin2 + 1):length(run2)]

post1 <- post1[seq(1, length(post1), by = thin_every)]
post2 <- post2[seq(1, length(post2), by = thin_every)]

trees <- c(post1, post2)
runs <- c(rep("run1", length(post1)), rep("run2", length(post2)))

message("After burn-in (and thinning every ", thin_every, "): N=", length(trees), " trees")

for (i in seq_along(trees)) {
  trees[[i]]$tip.label <- normalize_labels(trees[[i]]$tip.label)
}

# =========================
# 3) RF DISTANCES
# =========================

message("Computing RF distance matrix (may take time with many trees)...")
RF <- RF.dist(trees)
rf_dist <- as.dist(RF)
write.csv(as.matrix(RF), file.path(outdir, "RF_matrix_all_trees.csv"))

# =========================
# 4) CHOOSE OPTIMAL K (SILHOUETTE ON HIERARCHICAL CLUSTERING)
# =========================

# Hierarchical clustering on RF distances (Ward.D2)
hc <- hclust(rf_dist, method = "ward.D2")

message("Evaluating silhouette scores for K=", K_min, "..", K_max, " (cutree on Ward.D2 hclust)...")
Ks <- K_min:K_max

sil_scores <- sapply(Ks, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, rf_dist)[, 3])
})

sil_df <- data.frame(k = Ks, score = sil_scores)
K_opt <- sil_df$k[which.max(sil_df$score)]
message("Chosen K (max silhouette) = ", K_opt)

p_sil <- ggplot(sil_df, aes(k, score)) +
  geom_line() +
  geom_point(size = 2.5) +
  geom_vline(xintercept = K_opt, linetype = "dashed") +
  labs(
    title = "Silhouette scores for K (Ward.D2 hierarchical clustering on RF distances)",
    x = "Number of clusters (K)",
    y = "Mean silhouette width"
  ) +
  theme_report()

save_plot(p_sil, "01_silhouette_scores_vs_K.pdf")

# Final cluster assignment at optimal K
clusters <- factor(cutree(hc, k = K_opt)) # 1..K_opt
cluster_cols <- make_cluster_palette(K_opt)

write.csv(
  data.frame(tree_index = seq_along(trees), run = runs, cluster = clusters),
  file.path(outdir, "cluster_assignments_hclust_optK.csv"),
  row.names = FALSE
)

# =========================
# 5) MDS (2D) TREE SPACE
# =========================

mds2 <- cmdscale(rf_dist, k = 2)

mds_df <- data.frame(
  x = mds2[, 1],
  y = mds2[, 2],
  run = factor(runs, levels = c("run1", "run2")),
  cluster = clusters
)

# ---- Re-label clusters so they follow the visible left->right structure on MDS1 ----
# PAM cluster IDs are arbitrary (e.g. "1" and "2" may swap). To make plots consistent,
# we reassign cluster numbers by the mean MDS1 position: cluster 1 = leftmost centroid.
mds_df$cluster_raw <- as.integer(mds_df$cluster)
centroids <- aggregate(x ~ cluster_raw, data = mds_df, FUN = mean)
order_lr <- centroids$cluster_raw[order(centroids$x)] # left -> right on MDS1
map_lr <- setNames(seq_along(order_lr), order_lr)

# Update BOTH the global `clusters` and the column in `mds_df` so downstream steps match.
clusters <- factor(as.integer(map_lr[as.character(as.integer(clusters))]),
  levels = seq_len(K_opt)
)
mds_df$cluster <- clusters

# Overwrite cluster assignment table using the left->right relabeled clusters
write.csv(
  data.frame(tree_index = seq_along(trees), run = runs, cluster = clusters),
  file.path(outdir, "cluster_assignments_hclust_optK.csv"),
  row.names = FALSE
)

# 5a) MDS colored by run
p_mds_runs <- ggplot(mds_df, aes(x, y, color = run)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = run_cols, name = "Run") +
  labs(title = "Tree space (MDS on RF distances): colored by run", x = "MDS 1", y = "MDS 2") +
  theme_report()

save_plot(p_mds_runs, "02_MDS_colored_by_run.pdf")

# 5b) MDS colored by cluster (optimal K)
p_mds_cluster <- ggplot(mds_df, aes(x, y, color = cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = cluster_cols, name = paste0("Cluster (K=", K_opt, ")")) +
  labs(
    title = paste0("Tree space (MDS on RF distances): colored by cluster (K=", K_opt, ")"),
    x = "MDS 1", y = "MDS 2"
  ) +
  theme_report()

save_plot(p_mds_cluster, "03_MDS_colored_by_cluster.pdf")

# =========================
# 6) DENDROGRAM (Ward.D2) COLORED BY CLUSTER
# =========================

# Re-use the Ward.D2 clustering and color branches by K
dend <- as.dendrogram(hc)
labels(dend) <- rep("", length(labels(dend)))

dend_col <- dendextend::color_branches(dend, k = K_opt, col = unname(cluster_cols))

pdf(file.path(outdir, "04_dendrogram_branches_colored_by_cluster.pdf"),
  width = plot_width, height = plot_height
)
par(family = base_family, cex = 0.7, mar = c(5, 4, 4, 2) + 0.1)
plot(dend_col,
  main = paste0("Hierarchical clustering (Ward.D2): branches colored (K=", K_opt, ")"),
  ylab = "Height", xlab = ""
)
legend("topright",
  legend = paste0("Cluster ", seq_len(K_opt)),
  fill = unname(cluster_cols[as.character(seq_len(K_opt))]),
  border = NA, bty = "n", cex = 0.8
)
dev.off()

# =========================
# 7) CONSENSUS TREES PER (OPTIMAL) CLUSTER
# =========================

message("Writing consensus trees per cluster (p=", consensus_p, ")...")
cluster_cons <- vector("list", length = nlevels(clusters))

for (i in seq_len(nlevels(clusters))) {
  idx <- which(as.integer(clusters) == i)
  trees_i <- trees[idx]
  message("  Cluster ", i, ": ", length(trees_i), " trees")
  cons_i <- consensus(trees_i, p = consensus_p)
  cluster_cons[[i]] <- cons_i
  write.tree(cons_i, file.path(outdir, paste0("consensus_cluster_", i, "_K", K_opt, ".tre")))
}

# =========================
# 8) ML VS CLUSTER-CONSENSUS RF COMPARISONS
# =========================

message("Running ML vs cluster-consensus RF comparisons...")

ML <- read.tree(ml_tree_path)
ML <- drop_fake_numeric_tips(ML)
ML$tip.label <- normalize_labels(ML$tip.label)

# Normalize consensus labels
for (i in seq_along(cluster_cons)) {
  cluster_cons[[i]]$tip.label <- normalize_labels(cluster_cons[[i]]$tip.label)
}

# Keep only taxa common to ML and all cluster consensuses
common_tips <- Reduce(intersect, c(list(ML$tip.label), lapply(cluster_cons, function(t) t$tip.label)))
message("Common taxa across ML + all consensuses: ", length(common_tips))

# --- prune to common tips
ML2 <- keep_common_tips(ML, common_tips)
cons2 <- lapply(cluster_cons, \(t) keep_common_tips(t, common_tips))

# --- drop any NULLs (important!)
cons2 <- Filter(Negate(is.null), cons2)

# --- sanity checks (helps pinpoint which tree is bad)
stopifnot(inherits(ML2, "phylo"))
bad <- which(!vapply(cons2, inherits, logical(1), "phylo"))
if (length(bad)) {
  stop("These consensus trees are not 'phylo' after pruning: ", paste(bad, collapse = ", "))
}
nTips <- c(ML = ape::Ntip(ML2), vapply(cons2, ape::Ntip, integer(1)))
if (any(nTips < 3)) {
  stop("Some trees have < 3 tips after pruning: ", paste(names(nTips)[nTips < 3], collapse = ", "))
}

# --- build a proper multiPhylo (portable way)
comp_trees <- c(list(ML2), cons2)
class(comp_trees) <- "multiPhylo"

comp_names <- c("ML", paste0("cl", seq_along(cons2)))

rf_matrix <- as.matrix(phangorn::RF.dist(comp_trees))
dimnames(rf_matrix) <- list(comp_names, comp_names)

# Also output a simple "ML vs each cluster consensus" table
ml_vs_clusters <- data.frame(
  cluster = paste0("cl", seq_along(cons2)),
  RF_to_ML = rf_matrix["ML", paste0("cl", seq_along(cons2))]
)
write.csv(ml_vs_clusters, file.path(outdir, "RF_ML_vs_each_clusterConsensus.csv"), row.names = FALSE)

# ggplot heatmap (same theme/fonts/dimensions)
rf_long <- as.data.frame(rf_matrix) |>
  tibble::rownames_to_column("tree1") |>
  pivot_longer(-tree1, names_to = "tree2", values_to = "RF")

p_rf_heat <- ggplot(rf_long, aes(tree2, tree1, fill = RF)) +
  geom_tile() +
  coord_equal() +
  labs(
    title = "RF distance heatmap: ML and cluster consensus trees",
    x = NULL, y = NULL
  ) +
  theme_report() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_rf_heat, "05_RF_heatmap_ML_vs_clusterConsensus.pdf", w = plot_width, h = plot_height)

message("Done. Outputs written to: ", normalizePath(outdir))
