devtools::load_all()
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

ppi_edges_list <- list(
    all = ppi_edges, 
    high_conf = ppi_edges_high_conf)

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

for (i in names(graphs)) {
    saveRDS(graphs[[i]], file.path(out_dir, paste0(i, "_graph.rds")))
}

degree_graph <- graphs$all
weighted_degree_graph <- graphs$all
degree_high_conf_graph <- graphs$high_conf
weighted_degree_high_conf_graph <- graphs$high_conf

k <- 128

landmarks <- list(
  degree_all = gsemb_select_landmarks(degree_graph, k = k, method = "degree"),
  weighted_degree_all = gsemb_select_landmarks(weighted_degree_graph, k = k, method = "degree"),
  degree_high = gsemb_select_landmarks(degree_high_conf_graph, k = k, method = "degree"),
  weighted_degree_high = gsemb_select_landmarks(weighted_degree_high_conf_graph, k = k, method = "degree")
)

# -----------------------------------------------------------------------------------------------------------
# Calculate the landmark jaccard similarity
#两者值越大证明相似度越高
source("scripts/landmark_jaccard.R")
out1 <- landmark_jaccard(landmarks)
print(out1)
# dir.create(file.path(out_dir, "degree"), showWarnings = FALSE)
write.csv(out1, file.path(out_dir, "degree/degree_landmark_jaccard.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark spearman correlation （-1，1）
# 值越大证明两者越相关
scores <- list(
  degree_all = setNames(as.numeric(Matrix::rowSums(degree_graph != 0)), rownames(degree_graph)),
  weighted_degree_all = setNames(as.numeric(Matrix::rowSums(weighted_degree_graph)), rownames(weighted_degree_graph)),
  degree_high = setNames(as.numeric(Matrix::rowSums(degree_high_conf_graph != 0)), rownames(degree_high_conf_graph)),
  weighted_degree_high = setNames(as.numeric(Matrix::rowSums(weighted_degree_high_conf_graph)), rownames(weighted_degree_high_conf_graph))
)

source("scripts/landmark_spearman.R")
out2 <- landmark_spearman(scores)
print(out2)
write.csv(out2, file.path(out_dir, "degree/degree_landmark_spearman.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark quality
# 置信度越高 连接的节点越多 证明质量越高
source("scripts/landmark_quality.R")
out3 <- landmark_quality(landmarks, graphs$all, graphs$high_conf)
print(out3)
write.csv(out3, file.path(out_dir, "degree/degree_landmark_quality.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark coverage@2 and mean nearest landmark distance
# 覆盖率越高 证明质量越高 距离越短 证明质量越高
source("scripts/landmark_coverage.R")
out4 <- landmark_coverage(landmarks, graphs$high_conf, radius = 2L)
print(out4)
write.csv(out4, file.path(out_dir, "degree/degree_landmark_coverage.csv"))
# -----------------------------------------------------------------------------------------------------------
# calculate the landmark dropout 5% edge-drop Jaccard stability
# 值越大证明稳定性越高
source("scripts/landmark_dropout.R")
out5 <- landmark_dropout(landmarks, edges = ppi_edges, k = k,
                            dropout_fraction = 0.05, repeats = 100L, seed = 20260917L)
print(out5)
write.csv(out5, file.path(out_dir, "degree/degree_landmark_dropout.csv"))
# -----------------------------------------------------------------------------------------------------------
result <- "E:/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks/test_results"
