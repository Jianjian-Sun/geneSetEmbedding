devtools::load_all()
# -----------------------------------------------------------------------------------------------------------
# K-Medoids landmark geometry (farthest-first + restricted PAM)
# Four methods: kmedoids / weighted_kmedoids × all / high
# Skip: Spearman (no global node score), dropout, cross-Jaccard (later)
# -----------------------------------------------------------------------------------------------------------
ppi_links_file <- "E:/genesetembedding/step1/STRING/9606.protein.links.v12.0.txt"
out_dir <- "E:/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks/test_results"
kmedoids_dir <- file.path(out_dir, "kmedoids")
dir.create(kmedoids_dir, recursive = TRUE, showWarnings = FALSE)

k <- 128L
seed <- 20260919L
kmedoids_m <- 50L
kmedoids_max_iter <- 10L

all_graph_rds <- file.path(out_dir, "all_graph.rds")
high_graph_rds <- file.path(out_dir, "high_conf_graph.rds")

if (file.exists(all_graph_rds) && file.exists(high_graph_rds)) {
  graphs <- list(
    all = readRDS(all_graph_rds),
    high_conf = readRDS(high_graph_rds)
  )
} else {
  ppi_edges <- read.table(
    ppi_links_file,
    header = TRUE,
    sep = " ",
    stringsAsFactors = FALSE
  )
  ppi_edges_high_conf <- subset(ppi_edges, combined_score > 700)
  graphs <- lapply(
    list(all = ppi_edges, high_conf = ppi_edges_high_conf),
    function(edges) {
      gsemb_build_graph(
        edges,
        node1 = "protein1",
        node2 = "protein2",
        weight = "combined_score",
        directed = FALSE
      )
    }
  )
  for (i in names(graphs)) {
    saveRDS(graphs[[i]], file.path(out_dir, paste0(i, "_graph.rds")))
  }
}

select_km <- function(adj, method) {
  gsemb_select_landmarks(
    adj,
    k = k,
    method = method,
    seed = seed,
    kmedoids_m = kmedoids_m,
    kmedoids_max_iter = kmedoids_max_iter
  )
}

# -----------------------------------------------------------------------------------------------------------
# Select landmarks
# -----------------------------------------------------------------------------------------------------------
landmarks <- list(
  kmedoids_all = select_km(graphs$all, "kmedoids"),
  weighted_kmedoids_all = select_km(graphs$all, "kmedoids"),
  kmedoids_high = select_km(graphs$high_conf, "kmedoids"),
  weighted_kmedoids_high = select_km(graphs$high_conf, "kmedoids")
)

landmarks_meta <- list(
  k = k,
  seed = seed,
  kmedoids_m = kmedoids_m,
  kmedoids_max_iter = kmedoids_max_iter,
  created = Sys.time(),
  methods = names(landmarks)
)
saveRDS(
  list(landmarks = landmarks, meta = landmarks_meta),
  file.path(kmedoids_dir, "kmedoids_landmarks.rds")
)
landmarks_df <- as.data.frame(
  lapply(landmarks, function(x) {
    `length<-`(as.character(x), k)
  }),
  stringsAsFactors = FALSE
)
write.csv(
  landmarks_df,
  file.path(kmedoids_dir, "kmedoids_landmarks.csv"),
  row.names = FALSE
)
# Later:
#   obj <- readRDS(file.path(kmedoids_dir, "kmedoids_landmarks.rds"))
#   landmarks <- obj$landmarks

# -----------------------------------------------------------------------------------------------------------
# Jaccard (set overlap; higher = more similar)
# -----------------------------------------------------------------------------------------------------------
source("scripts/landmark_jaccard.R")
out1 <- landmark_jaccard(landmarks)
print(out1)
write.csv(out1, file.path(kmedoids_dir, "kmedoids_landmark_jaccard.csv"))

# -----------------------------------------------------------------------------------------------------------
# Spearman: skipped — K-Medoids has no global per-node score like degree/betweenness
# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------
# Quality (hub / high-confidence character of selected landmarks)
# -----------------------------------------------------------------------------------------------------------
source("scripts/landmark_quality.R")
out3 <- landmark_quality(landmarks, graphs$all, graphs$high_conf)
print(out3)
write.csv(out3, file.path(kmedoids_dir, "kmedoids_landmark_quality.csv"))

# -----------------------------------------------------------------------------------------------------------
# Coverage@2 on high graph (higher coverage / lower mean distance = better spread)
# -----------------------------------------------------------------------------------------------------------
source("scripts/landmark_coverage.R")
out4 <- landmark_coverage(landmarks, graphs$high_conf, radius = 2L)
print(out4)
write.csv(out4, file.path(kmedoids_dir, "kmedoids_landmark_coverage.csv"))

# -----------------------------------------------------------------------------------------------------------
# Dropout / cross-Jaccard: later
# -----------------------------------------------------------------------------------------------------------
