# Helpers for the SLE Landmark real-cohort benchmark.
#
# This file intentionally contains only glue that is not already provided by
# geneSetEmbedding or step3/step3_modular. It does not reimplement RWR, SVD,
# Gaussian fitting, Weighted-GSEA, ORA, Jaccard, or Spearman correlation.

load_string_symbol_map <- function(path) {
  info <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("string_protein_id", "preferred_name")
  missing_cols <- setdiff(required, names(info))
  if (length(missing_cols)) {
    stop("STRING info is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  info$string_protein_id <- as.character(info$string_protein_id)
  info$preferred_name <- as.character(info$preferred_name)
  keep <- !is.na(info$string_protein_id) & nzchar(info$string_protein_id) &
    !is.na(info$preferred_name) & nzchar(info$preferred_name)
  info <- info[keep, , drop = FALSE]
  if (anyDuplicated(info$string_protein_id)) {
    stop("STRING info has duplicated string_protein_id values.", call. = FALSE)
  }
  info
}

map_ensp_to_symbols <- function(ids, protein_info) {
  ids <- as.character(ids)
  m <- match(ids, protein_info$string_protein_id)
  if (anyNA(m)) {
    stop(
      "Unmapped STRING IDs: ",
      paste(head(ids[is.na(m)], 10L), collapse = ", "),
      call. = FALSE
    )
  }
  protein_info$preferred_name[m]
}

remap_adj_to_symbols <- function(adj, protein_info) {
  ids <- rownames(adj)
  if (is.null(ids) || is.null(colnames(adj))) {
    stop("adj must have rownames and colnames.", call. = FALSE)
  }
  if (!identical(ids, colnames(adj))) {
    stop("adj rownames and colnames must match.", call. = FALSE)
  }
  sym <- map_ensp_to_symbols(ids, protein_info)
  if (anyDuplicated(sym)) {
    # ponytail: STRING v12 preferred_name is unique for this graph; collapse only if needed
    strength <- as.numeric(Matrix::rowSums(adj))
    ord <- order(strength, decreasing = TRUE)
    keep <- sort(ord[!duplicated(sym[ord])])
    adj <- adj[keep, keep, drop = FALSE]
    sym <- sym[keep]
  }
  dimnames(adj) <- list(sym, sym)
  adj
}

prepare_ranked_stats <- function(deg_df,
                                 gene_col = "gene_symbol",
                                 stat_col = "log2FoldChange") {
  missing_cols <- setdiff(c(gene_col, stat_col), names(deg_df))
  if (length(missing_cols)) {
    stop("DEG table is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  genes <- as.character(deg_df[[gene_col]])
  values <- as.numeric(deg_df[[stat_col]])
  keep <- !is.na(genes) & nzchar(genes) & is.finite(values)
  ranked <- data.frame(gene = genes[keep], value = values[keep], stringsAsFactors = FALSE)
  if (!nrow(ranked)) stop("DEG table has no usable ranked statistics.", call. = FALSE)
  ranked <- ranked[order(ranked$gene, -abs(ranked$value)), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$gene), , drop = FALSE]
  out <- ranked$value
  names(out) <- ranked$gene
  out
}

build_fit_from_landmarks <- function(adj,
                                     gene_sets,
                                     landmarks,
                                     dim = 64L,
                                     alpha = 0.5,
                                     tol = 1e-10,
                                     max_iter = 200L,
                                     normalize = "col",
                                     seed = 1L) {
  nodes <- rownames(adj)
  if (is.null(nodes)) stop("adj must have node names.", call. = FALSE)
  landmarks <- unique(as.character(landmarks))
  missing_landmarks <- setdiff(landmarks, nodes)
  if (length(missing_landmarks)) {
    stop("Landmarks absent from graph: ", paste(head(missing_landmarks, 10L), collapse = ", "), call. = FALSE)
  }
  if (!length(landmarks)) stop("landmarks must be non-empty.", call. = FALSE)

  node_features <- gsemb_compute_node_landmark_features(
    adj = adj,
    landmarks = landmarks,
    alpha = alpha,
    tol = tol,
    max_iter = max_iter,
    normalize = normalize,
    seed = seed
  )
  max_dim <- min(nrow(node_features), ncol(node_features))
  if (dim < 1L || dim > max_dim) {
    stop("dim must be between 1 and ", max_dim, " for these features.", call. = FALSE)
  }
  emb <- gsemb_fit_gene_embedding(
    node_features = node_features,
    dim = as.integer(dim),
    method = "svd",
    seed = seed
  )$embedding
  gauss <- gsemb_fit_set_gaussians_from_members(emb, gene_sets)
  fit <- list(
    adj = adj,
    method = "svd",
    gene_embedding = emb,
    set_mu = gauss$mu,
    set_var = gauss$var,
    landmarks = landmarks,
    losses = NULL
  )
  class(fit) <- "gsemb_embedding"
  fit
}

as_comparable_enrichment <- function(x) {
  required <- c("ID", "pvalue", "p.adjust")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols)) {
    stop("Enrichment result is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if ("status" %in% names(x)) {
    x <- x[!is.na(x$status) & x$status == "ok", , drop = FALSE]
  }
  x <- x[!is.na(x$ID) & nzchar(as.character(x$ID)), , drop = FALSE]
  x[order(x$p.adjust, x$pvalue, na.last = TRUE), , drop = FALSE]
}

summarize_method_reproducibility <- function(pairwise_df, route) {
  required <- c("version", "pathway_mode", "jaccard", "overlap_coefficient", "rank_spearman")
  missing_cols <- setdiff(required, names(pairwise_df))
  if (length(missing_cols)) {
    stop("Pairwise table is missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  keys <- interaction(pairwise_df$version, pairwise_df$pathway_mode, drop = TRUE)
  rows <- lapply(split(pairwise_df, keys), function(d) {
    data.frame(
      route = route,
      method = as.character(d$version[[1]]),
      pathway_mode = as.character(d$pathway_mode[[1]]),
      n_cohort_pairs = nrow(d),
      mean_jaccard = mean(d$jaccard, na.rm = TRUE),
      mean_overlap_coefficient = mean(d$overlap_coefficient, na.rm = TRUE),
      mean_rank_spearman = mean(d$rank_spearman, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

load_landmark_candidates <- function(degree_rds, betweenness_rds, kmedoids_rds) {
  degree <- readRDS(degree_rds)$landmarks
  betweenness <- readRDS(betweenness_rds)$landmarks
  kmedoids <- readRDS(kmedoids_rds)$landmarks
  out <- list(
    weighted_degree_high = degree$weighted_degree_high,
    weighted_betweenness_high = betweenness$weighted_betweenness_high,
    weighted_kmedoids_high = kmedoids$weighted_kmedoids_high
  )
  missing <- names(out)[vapply(out, is.null, logical(1))]
  if (length(missing)) {
    stop("Missing selected Landmark sets: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out
}

write_run_metadata <- function(path, values) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- paste0(names(values), ": ", vapply(values, paste, collapse = ", ", FUN.VALUE = character(1)))
  writeLines(lines, path)
  invisible(path)
}
