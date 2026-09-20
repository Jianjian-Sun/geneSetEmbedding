devtools::load_all()
# -----------------------------------------------------------------------------------------------------------
# Load the PPI edges and build the graphs/select landmarks
ppi_links_file <- "E:/genesetembedding/step1/STRING/9606.protein.links.v12.0.txt"
ppi_edges <- read.table(ppi_links_file, header = TRUE, sep = " ", stringsAsFactors = FALSE)
head(ppi_edges)
ppi_edges_high_conf <- subset(ppi_edges, combined_score > 700)

ppi_edges_list <- list(
  all = ppi_edges,
  high_conf = ppi_edges_high_conf
)

graphs <- lapply(ppi_edges_list, function(edges) {
  gsemb_build_graph(
    edges,
    node1 = "protein1",
    node2 = "protein2",
    weight = "combined_score",
    directed = FALSE
  )
})

out_dir <- "E:/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks/test_results"
betweenness_dir <- file.path(out_dir, "betweenness")
dir.create(betweenness_dir, recursive = TRUE, showWarnings = FALSE)

for (i in names(graphs)) {
  saveRDS(graphs[[i]], file.path(out_dir, paste0(i, "_graph.rds")))
}

betweenness_graph <- graphs$all
weighted_betweenness_graph <- graphs$all
betweenness_high_conf_graph <- graphs$high_conf
weighted_betweenness_high_conf_graph <- graphs$high_conf

k <- 128L
# package default; kept explicit for the script
betweenness_cutoff <- 4

landmarks <- list(
  betweenness_all = gsemb_select_landmarks(
    betweenness_graph, k = k, method = "betweenness",
    betweenness_cutoff = betweenness_cutoff
  ),
  weighted_betweenness_all = gsemb_select_landmarks(
    weighted_betweenness_graph, k = k, method = "weighted_betweenness",
    betweenness_cutoff = betweenness_cutoff
  ),
  betweenness_high = gsemb_select_landmarks(
    betweenness_high_conf_graph, k = k, method = "betweenness",
    betweenness_cutoff = betweenness_cutoff
  ),
  weighted_betweenness_high = gsemb_select_landmarks(
    weighted_betweenness_high_conf_graph, k = k, method = "weighted_betweenness",
    betweenness_cutoff = betweenness_cutoff
  )
)
# -----------------------------------------------------------------------------------------------------------
# Persist landmarks for later reuse (skip re-running betweenness)
landmarks_meta <- list(
  k = k,
  betweenness_cutoff = betweenness_cutoff,
  created = Sys.time(),
  methods = names(landmarks)
)
saveRDS(
  list(landmarks = landmarks, meta = landmarks_meta),
  file.path(betweenness_dir, "betweenness_landmarks.rds")
)
# Optional plain-text dump (one column per method, unequal lengths padded with NA)
landmarks_df <- as.data.frame(
  lapply(landmarks, function(x) {
    `length<-`(as.character(x), k)
  }),
  stringsAsFactors = FALSE
)
write.csv(
  landmarks_df,
  file.path(betweenness_dir, "betweenness_landmarks.csv"),
  row.names = FALSE
)
# Later:
#   obj <- readRDS(".../test_results/betweenness/betweenness_landmarks.rds")
#   landmarks <- obj$landmarks
# -----------------------------------------------------------------------------------------------------------
# Full-node scores for Spearman (same formulas as gsemb_select_landmarks)
# 这一步暂时不做吧 介数中心算这个好像不划算 时间太长了
betweenness_scores <- function(adj, weighted = FALSE, cutoff = 4) {
  nodes <- rownames(adj)
  if (is.null(nodes)) stop("adj must have rownames")
  if (!weighted) {
    g <- igraph::graph_from_adjacency_matrix(
      adj != 0, mode = "undirected", diag = FALSE
    )
    bc <- igraph::betweenness(g, directed = FALSE, cutoff = cutoff)
  } else {
    g <- igraph::graph_from_adjacency_matrix(
      adj, mode = "undirected", weighted = TRUE, diag = FALSE
    )
    costs <- 1 / pmax(igraph::E(g)$weight, .Machine$double.eps)
    bc <- igraph::betweenness(
      g, directed = FALSE, weights = costs, cutoff = cutoff
    )
  }
  setNames(as.numeric(bc), nodes)
}

# -----------------------------------------------------------------------------------------------------------
# Calculate the landmark jaccard similarity
# 两者值越大证明相似度越高
source("scripts/landmark_jaccard.R")
out1 <- landmark_jaccard(landmarks)
print(out1)
write.csv(out1, file.path(betweenness_dir, "betweenness_landmark_jaccard.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark spearman correlation （-1，1）
# 值越大证明两者越相关
scores <- list(
  betweenness_all = betweenness_scores(
    betweenness_graph, weighted = FALSE, cutoff = betweenness_cutoff
  ),
  weighted_betweenness_all = betweenness_scores(
    weighted_betweenness_graph, weighted = TRUE, cutoff = betweenness_cutoff
  ),
  betweenness_high = betweenness_scores(
    betweenness_high_conf_graph, weighted = FALSE, cutoff = betweenness_cutoff
  ),
  weighted_betweenness_high = betweenness_scores(
    weighted_betweenness_high_conf_graph, weighted = TRUE, cutoff = betweenness_cutoff
  )
)

source("scripts/landmark_spearman.R")
out2 <- landmark_spearman(scores)
print(out2)
write.csv(out2, file.path(betweenness_dir, "betweenness_landmark_spearman.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark quality
# 置信度越高 连接的节点越多 证明质量越高
source("scripts/landmark_quality.R")
out3 <- landmark_quality(landmarks, graphs$all, graphs$high_conf)
print(out3)
write.csv(out3, file.path(betweenness_dir, "betweenness_landmark_quality.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark coverage@2 and mean nearest landmark distance
# 覆盖率越高 证明质量越高 距离越短 证明质量越高
source("scripts/landmark_coverage.R")
out4 <- landmark_coverage(landmarks, graphs$high_conf, radius = 2L)
print(out4)
write.csv(out4, file.path(betweenness_dir, "betweenness_landmark_coverage.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark dropout 5% edge-drop Jaccard stability
# 值越大证明稳定性越高
# repeats=20（非 degree 的 100）：每次需重算介数，全边图上 100 次过慢
source("scripts/landmark_dropout.R")
selectors <- list(
  betweenness_all = list(graph = "all", method = "betweenness"),
  weighted_betweenness_all = list(graph = "all", method = "weighted_betweenness"),
  betweenness_high = list(graph = "high", method = "betweenness"),
  weighted_betweenness_high = list(graph = "high", method = "weighted_betweenness")
)
out5 <- landmark_dropout(
  landmarks,
  edges = ppi_edges,
  k = k,
  dropout_fraction = 0.05,
  repeats = 20L,
  seed = 20260918L,
  selectors = selectors
)
print(out5$summary)
write.csv(out5$summary, file.path(betweenness_dir, "betweenness_landmark_dropout.csv"))
# -----------------------------------------------------------------------------------------------------------
