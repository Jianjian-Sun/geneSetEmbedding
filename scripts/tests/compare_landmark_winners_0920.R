# Pairwise Jaccard + 3-set Venn for family winners on high graph.
# Arms: weighted_degree_high, weighted_betweenness_high, weighted_kmedoids_high
devtools::load_all(".", quiet = TRUE)
source("scripts/landmark_jaccard.R")

root <- "E:/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks"
out_dir <- file.path(root, "test_results", "cross_method")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# degree run never saved landmarks RDS; recompute Top-K (cheap).
adj_high <- readRDS(file.path(root, "test_results", "high_conf_graph.rds"))
wd_high <- gsemb_select_landmarks(adj_high, k = 128L, method = "weighted_degree")

bet <- readRDS(file.path(root, "test_results", "betweenness", "betweenness_landmarks.rds"))
km <- readRDS(file.path(root, "test_results", "kmedoids", "kmedoids_landmarks.rds"))

landmarks <- list(
  weighted_degree_high = wd_high,
  weighted_betweenness_high = bet$landmarks$weighted_betweenness_high,
  weighted_kmedoids_high = km$landmarks$weighted_kmedoids_high
)
stopifnot(all(lengths(landmarks) == 128L))

jacc <- landmark_jaccard(landmarks)
write.csv(jacc, file.path(out_dir, "winners_pairwise_jaccard.csv"), row.names = FALSE)
print(jacc)

A <- landmarks$weighted_degree_high
B <- landmarks$weighted_betweenness_high
C <- landmarks$weighted_kmedoids_high
lab <- c("weighted_degree_high", "weighted_betweenness_high", "weighted_kmedoids_high")

only_a <- setdiff(A, union(B, C))
only_b <- setdiff(B, union(A, C))
only_c <- setdiff(C, union(A, B))
ab <- setdiff(intersect(A, B), C)
ac <- setdiff(intersect(A, C), B)
bc <- setdiff(intersect(B, C), A)
abc <- Reduce(intersect, list(A, B, C))

venn_counts <- data.frame(
  region = c("only_degree", "only_betweenness", "only_kmedoids",
             "degree_betweenness", "degree_kmedoids", "betweenness_kmedoids",
             "all_three"),
  n = c(length(only_a), length(only_b), length(only_c),
        length(ab), length(ac), length(bc), length(abc))
)
write.csv(venn_counts, file.path(out_dir, "winners_venn_counts.csv"), row.names = FALSE)
print(venn_counts)

# ponytail: base graphics 3-circle Venn; upgrade to ggVennDiagram if layout needs polish.
png(file.path(out_dir, "winners_venn.png"), width = 900, height = 800, res = 120)
par(mar = c(1, 1, 3, 1))
plot(NA, xlim = c(-1.6, 1.6), ylim = c(-1.5, 1.5), axes = FALSE, xlab = "", ylab = "",
     main = "Landmark winners (k=128, high graph)")
theta <- seq(0, 2 * pi, length.out = 200)
cx <- c(-0.55, 0.55, 0)
cy <- c(0.25, 0.25, -0.55)
cols <- c("#4C78A8", "#F58518", "#54A24B")
for (i in 1:3) {
  polygon(cx[i] + 1.05 * cos(theta), cy[i] + 1.05 * sin(theta),
          border = cols[i], lwd = 2, col = adjustcolor(cols[i], 0.18))
}
text(-1.15, 1.15, lab[1], col = cols[1], cex = 0.85, font = 2)
text(1.15, 1.15, lab[2], col = cols[2], cex = 0.85, font = 2)
text(0, -1.35, lab[3], col = cols[3], cex = 0.85, font = 2)
text(-0.85, 0.35, length(only_a), cex = 1.1, font = 2)
text(0.85, 0.35, length(only_b), cex = 1.1, font = 2)
text(0, -0.95, length(only_c), cex = 1.1, font = 2)
text(0, 0.45, length(ab), cex = 1.05, font = 2)
text(-0.45, -0.25, length(ac), cex = 1.05, font = 2)
text(0.45, -0.25, length(bc), cex = 1.05, font = 2)
text(0, -0.05, length(abc), cex = 1.2, font = 2)
dev.off()

saveRDS(
  list(landmarks = landmarks, jaccard = jacc, venn_counts = venn_counts,
       created = Sys.time()),
  file.path(out_dir, "winners_landmarks.rds")
)

cat("Wrote:\n",
    file.path(out_dir, "winners_pairwise_jaccard.csv"), "\n",
    file.path(out_dir, "winners_venn_counts.csv"), "\n",
    file.path(out_dir, "winners_venn.png"), "\n", sep = "")
