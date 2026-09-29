# Time weighted landmark selection on STRING combined_score > 700.
# Methods: degree, betweenness, kmedoids.
# Runs are sequential so each elapsed time is that method alone.
#
# From the weighted-degree-landmarks worktree root:
#   Rscript scripts/tests/time_weighted_landmarks_score700.R


# debug = FALSE compiles the C++ at -O2. devtools::load_all() drops this argument.
pkgload::load_all(debug = FALSE, recompile = TRUE)
# -----------------------------------------------------------------------------------------------------------
# Load the PPI edges and build the graphs/select landmarks
ppi_links_file <- "/mnt/d/genesetembedding/step1/STRING/9606.protein.links.v12.0.txt"
ppi_edges <- read.table(ppi_links_file, header = TRUE, sep = " ", stringsAsFactors = FALSE)
#head(ppi_edges)
ppi_edges <- ppi_edges[ppi_edges$protein1 != ppi_edges$protein2, ] # 删除自环 得到的结果是0个自环
ppi_edges <- ppi_edges[ppi_edges$protein1 < ppi_edges$protein2, ] #去掉反向重复边 保留了6857702条边
#nrow(ppi_edges)
ppi_edges_high_conf <- subset(ppi_edges, combined_score > 700) # 过滤低置信度边 472000条边(没去掉的时候) -> 236000条边(去掉反向重复边之后)
#nrow(ppi_edges_high_conf)

graphs <- gsemb_build_graph(
    ppi_edges_high_conf,
    node1 = "protein1",
    node2 = "protein2",
    weight = "combined_score",
    directed = FALSE
)
data_dir <- "/mnt/d/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks/test_results/data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(graphs, file.path(data_dir, "high_conf_graph_rectify0929.rds"))
# -----------------------------------------------------------------------------------------------------------

score_cutoff <- 700
k <- 128
seed <- 20260929
betweenness_cutoff <- -1
kmedoids_m <- 50
kmedoids_max_iter <- 10

methods <- c("degree", "betweenness", "kmedoids")
rows <- vector("list", length(methods))
landmark_rows <- vector("list", length(methods))

for (i in seq_along(methods)) {
  method <- methods[[i]]
  message("timing ", method)
  landmarks <- NULL
  elapsed <- system.time({
    landmarks <- gsemb_select_landmarks(
      graphs,
      k = k,
      method = method,
      seed = seed,
      betweenness_cutoff = betweenness_cutoff,
      kmedoids_m = kmedoids_m,
      kmedoids_max_iter = kmedoids_max_iter
    )
  })
  rows[[i]] <- data.frame(
    method = method,
    score_cutoff = score_cutoff,
    n_nodes = nrow(graphs),
    n_edges = Matrix::nnzero(graphs) / 2,
    k = length(landmarks),
    user_sec = unname(elapsed[["user.self"]]),
    system_sec = unname(elapsed[["sys.self"]]),
    elapsed_sec = unname(elapsed[["elapsed"]]),
    stringsAsFactors = FALSE
  )
  landmark_rows[[i]] <- data.frame(
    method = method,
    rank = seq_along(landmarks),
    landmark = landmarks,
    stringsAsFactors = FALSE
  )
  message(sprintf("  %s elapsed %.1f sec", method, elapsed[["elapsed"]]))
}


timing_dir <- "/mnt/d/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks/test_results/timing"
timing <- do.call(rbind, rows)
landmarks_out <- do.call(rbind, landmark_rows)
write.csv(timing, file.path(timing_dir, "weighted_landmarks_score700_timing_rectify0929.csv"), row.names = FALSE)
write.csv(landmarks_out, file.path(timing_dir, "landmarks_top128_rectify0929.csv"), row.names = FALSE)
