devtools::load_all()
# -----------------------------------------------------------------------------------------------------------
# Sensitivity: combined_score > 600 (mid) vs > 700 (high)
# Question: does mid separate degree vs weighted_degree more than high (Jaccard 0.94)?
# Skip: dropout (slow; add when mid Jaccard looks interesting)
# -----------------------------------------------------------------------------------------------------------
ppi_links_file <- "E:/genesetembedding/step1/STRING/9606.protein.links.v12.0.txt"
out_dir <- "E:/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks/test_results"
mid_dir <- file.path(out_dir, "degree_mid600")
dir.create(mid_dir, recursive = TRUE, showWarnings = FALSE)

k <- 128L
mid_thr <- 600
high_thr <- 700

ppi_edges <- read.table(
  ppi_links_file,
  header = TRUE,
  sep = " ",
  stringsAsFactors = FALSE
)
ppi_edges_mid <- subset(ppi_edges, combined_score > mid_thr)
ppi_edges_high <- subset(ppi_edges, combined_score > high_thr)

message(sprintf(
  "edges: all=%d mid(>%d)=%d high(>%d)=%d",
  nrow(ppi_edges), mid_thr, nrow(ppi_edges_mid), high_thr, nrow(ppi_edges_high)
))

build_g <- function(edges) {
  gsemb_build_graph(
    edges,
    node1 = "protein1",
    node2 = "protein2",
    weight = "combined_score",
    directed = FALSE
  )
}

high_rds <- file.path(out_dir, "high_conf_graph.rds")
graphs <- list(
  mid = build_g(ppi_edges_mid),
  high = if (file.exists(high_rds)) readRDS(high_rds) else build_g(ppi_edges_high)
)
saveRDS(graphs$mid, file.path(mid_dir, "mid600_graph.rds"))

landmarks <- list(
  degree_mid = gsemb_select_landmarks(graphs$mid, k = k, method = "degree"),
  weighted_degree_mid = gsemb_select_landmarks(graphs$mid, k = k, method = "degree"),
  degree_high = gsemb_select_landmarks(graphs$high, k = k, method = "degree"),
  weighted_degree_high = gsemb_select_landmarks(graphs$high, k = k, method = "degree")
)

saveRDS(
  list(
    landmarks = landmarks,
    meta = list(k = k, mid_thr = mid_thr, high_thr = high_thr, created = Sys.time())
  ),
  file.path(mid_dir, "degree_mid600_landmarks.rds")
)
write.csv(
  as.data.frame(lapply(landmarks, function(x) `length<-`(as.character(x), k)),
                stringsAsFactors = FALSE),
  file.path(mid_dir, "degree_mid600_landmarks.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------------------------------------
# Jaccard (focus: degree_mid vs weighted_degree_mid; high pair is the known ~0.94 anchor)
# -----------------------------------------------------------------------------------------------------------
source("scripts/landmark_jaccard.R")
out1 <- landmark_jaccard(landmarks)
print(out1)
write.csv(out1, file.path(mid_dir, "degree_mid600_landmark_jaccard.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------------------------------------
# Spearman on full-node scores (within-graph only; cross-graph scores not aligned)
# -----------------------------------------------------------------------------------------------------------
scores <- list(
  degree_mid = setNames(as.numeric(Matrix::rowSums(graphs$mid != 0)), rownames(graphs$mid)),
  weighted_degree_mid = setNames(as.numeric(Matrix::rowSums(graphs$mid)), rownames(graphs$mid)),
  degree_high = setNames(as.numeric(Matrix::rowSums(graphs$high != 0)), rownames(graphs$high)),
  weighted_degree_high = setNames(as.numeric(Matrix::rowSums(graphs$high)), rownames(graphs$high))
)
source("scripts/landmark_spearman.R")
out2 <- landmark_spearman(scores)
print(out2)
write.csv(out2, file.path(mid_dir, "degree_mid600_landmark_spearman.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------------------------------------
# Quality / Coverage@2 on high graph (same eval graph as degree_0917 for comparability)
# -----------------------------------------------------------------------------------------------------------
all_rds <- file.path(out_dir, "all_graph.rds")
adj_all <- if (file.exists(all_rds)) {
  readRDS(all_rds)
} else {
  build_g(ppi_edges)
}
source("scripts/landmark_quality.R")
out3 <- landmark_quality(landmarks, adj_all, graphs$high)
print(out3)
write.csv(out3, file.path(mid_dir, "degree_mid600_landmark_quality.csv"), row.names = FALSE)

source("scripts/landmark_coverage.R")
out4 <- landmark_coverage(landmarks, graphs$high, radius = 2L)
print(out4)
write.csv(out4, file.path(mid_dir, "degree_mid600_landmark_coverage.csv"), row.names = FALSE)

message("done -> ", mid_dir)
message("key row: degree_mid vs weighted_degree_mid (compare to degree_high vs weighted_degree_high ~0.94)")
