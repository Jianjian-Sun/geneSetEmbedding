#' Fit gene and gene set embeddings from PPI and gene sets
#'
#' This is the main entry point: build a graph from a PPI edge list (or accept
#' a pre-built adjacency matrix), compute diffusion features, learn a gene
#' embedding, and derive diagonal Gaussian embeddings for gene sets. Optionally,
#' train a Set2Gaussian model with torch.
#'
#' @param ppi Either a data.frame edge list or a Matrix adjacency matrix.
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param node1,node2 Column names in `ppi` when `ppi` is a data.frame.
#' @param weight Optional weight column name in `ppi` when `ppi` is a data.frame.
#' @param directed Logical; whether the graph is directed.
#' @param nodes Optional node ID vector used to restrict/order the graph.
#' @param method Embedding method: `"svd"`, `"torch_autoencoder"`, or `"set2gaussian_torch"`.
#' @param dim Embedding dimension.
#' @param k Number of landmarks used for diffusion features.
#' @param alpha Restart probability for diffusion.
#' @param tol Convergence tolerance for diffusion.
#' @param max_iter Maximum diffusion iterations.
#' @param normalize Normalization for `gsemb_transition_matrix`.
#' @param epochs,lr,batch_size Training hyperparameters for torch methods.
#' @param seed Random seed.
#' @param device `"cpu"` or `"cuda"` (falls back to CPU if CUDA is unavailable).
#' @param ... Passed through to the selected training backend.
#'
#' @return A `gsemb_embedding` object with components:
#'   `gene_embedding`, `set_mu`, `set_var`, `adj`, and training metadata.
#' @examples
#' \dontrun{
#' # Simulate a small PPI edge list
#' set.seed(42)
#' nodes <- paste0("GENE", 1:50)
#' edges <- data.frame(
#'   node1 = sample(nodes, 200, replace = TRUE),
#'   node2 = sample(nodes, 200, replace = TRUE),
#'   weight = runif(200, 0.5, 1)
#' )
#' # Remove self-loops
#' edges <- edges[edges$node1 != edges$node2, ]
#'
#' # Simulate a few gene sets
#' gene_sets <- list(
#'   SET1 = sample(nodes, 10),
#'   SET2 = sample(nodes, 8),
#'   SET3 = sample(nodes, 12)
#' )
#'
#' # Fit embedding using SVD (fastest method)
#' fit <- gsemb_fit(
#'   ppi = edges,
#'   gene_sets = gene_sets,
#'   method = "svd",
#'   dim = 16,
#'   k = 20,
#'   alpha = 0.5,
#'   epochs = 5
#' )
#'
#' # Inspect the result
#' print(fit)
#'
#' # Compute gene-gene similarities
#' sim <- gsemb_gene_similarity(fit)
#' dim(sim)
#' }
#' @export
gsemb_fit <- function(ppi,
                      gene_sets,
                      node1 = "node1",
                      node2 = "node2",
                      weight = NULL,
                      directed = FALSE,
                      nodes = NULL,
                      method = c("svd", "torch_autoencoder", "set2gaussian_torch"),
                      dim = 64,
                      k = 128,
                      alpha = 0.5,
                      tol = 1e-10,
                      max_iter = 200,
                      normalize = "col",
                      epochs = 200,
                      lr = 5e-3,
                      batch_size = 8,
                      seed = 1,
                      device = c("cpu", "cuda"),
                      ...) {
  method <- match.arg(method)
  device <- match.arg(device)
  gene_sets <- validate_gene_sets(gene_sets)

  adj <- if (inherits(ppi, "Matrix")) {
    ppi
  } else {
    gsemb_build_graph(ppi, node1 = node1, node2 = node2, weight = weight, directed = directed, nodes = nodes)
  }
  if (is.null(rownames(adj))) stop("PPI graph must have node names")

  if (method == "set2gaussian_torch") {
    fit <- gsemb_fit_set2gaussian_torch(
      adj = adj,
      gene_sets = gene_sets,
      k = k,
      dim = dim,
      alpha = alpha,
      tol = tol,
      max_iter = max_iter,
      normalize = normalize,
      epochs = epochs,
      lr = lr,
      batch_size = batch_size,
      seed = seed,
      device = device,
      ...
    )
    res <- list(
      adj = adj,
      method = method,
      gene_embedding = fit$gene_embedding,
      set_mu = fit$set_mu,
      set_var = fit$set_var,
      landmarks = fit$landmarks,
      losses = fit$losses
    )
    class(res) <- "gsemb_embedding"
    return(res)
  }

  landmarks <- gsemb_select_landmarks(adj, k = k, method = "degree", seed = seed)
  node_features <- gsemb_compute_node_landmark_features(
    adj = adj,
    landmarks = landmarks,
    alpha = alpha,
    tol = tol,
    max_iter = max_iter,
    normalize = normalize,
    seed = seed
  )
  emb_fit <- gsemb_fit_gene_embedding(
    node_features = node_features,
    dim = dim,
    method = if (method == "torch_autoencoder") "torch_autoencoder" else "svd",
    epochs = epochs,
    lr = lr,
    batch_size = max(64, 4 * dim),
    seed = seed,
    device = device
  )
  gene_embedding <- emb_fit$embedding
  gauss <- gsemb_fit_set_gaussians_from_members(gene_embedding, gene_sets)

  res <- list(
    adj = adj,
    method = method,
    gene_embedding = gene_embedding,
    set_mu = gauss$mu,
    set_var = gauss$var,
    landmarks = landmarks,
    losses = NULL
  )
  class(res) <- "gsemb_embedding"
  res
}

#' Compute gene-gene cosine similarities
#'
#' @param x A `gsemb_embedding` object or a numeric gene embedding matrix.
#' @param genes Optional character vector of genes to include as rows.
#' @param other_genes Optional character vector of genes to include as columns.
#' @param eps Small constant for numerical stability.
#'
#' @return A numeric similarity matrix.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' sim_all <- gsemb_gene_similarity(fit)
#'
#' # Compute similarities for a subset of genes
#' sim_sub <- gsemb_gene_similarity(fit, genes = c("GENE1", "GENE2", "GENE3"))
#' }
#' @export
gsemb_gene_similarity <- function(x,
                                  genes = NULL,
                                  other_genes = NULL,
                                  eps = 1e-12) {
  gene_emb <- if (inherits(x, "gsemb_embedding")) x$gene_embedding else x
  if (is.null(genes) && is.null(other_genes)) {
    return(gsemb_gene_cosine_similarity(gene_emb, eps = eps))
  }
  all_genes <- rownames(gene_emb)
  if (is.null(all_genes)) stop("gene embedding must have rownames")
  if (is.null(genes)) genes <- all_genes
  genes <- intersect(as.character(genes), all_genes)
  if (length(genes) == 0) stop("no genes found in embedding")
  if (is.null(other_genes)) {
    other_genes <- genes
  } else {
    other_genes <- intersect(as.character(other_genes), all_genes)
    if (length(other_genes) == 0) stop("no other_genes found in embedding")
  }
  gsemb_gene_cosine_similarity(gene_emb[genes, , drop = FALSE], gene_emb[other_genes, , drop = FALSE], eps = eps)
}

#' Compute set-set distances between Gaussian embeddings
#'
#' @param x A `gsemb_embedding` object.
#' @param sets Optional character vector of sets to include as rows.
#' @param other_sets Optional character vector of sets to include as columns.
#' @param metric Distance metric: `"w2"` or `"sym_kl"`.
#' @param eps Lower bound for variances.
#'
#' @return A numeric distance matrix.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Compute Wasserstein-2 distances between all gene sets
#' dist_w2 <- gsemb_set_similarity(fit, metric = "w2")
#'
#' # Compute symmetric KL divergences for a subset of sets
#' dist_kl <- gsemb_set_similarity(
#'   fit,
#'   sets = c("SET1", "SET2"),
#'   metric = "sym_kl"
#' )
#' }
#' @export
gsemb_set_similarity <- function(x,
                                 sets = NULL,
                                 other_sets = NULL,
                                 metric = c("w2", "sym_kl"),
                                 eps = 1e-8) {
  metric <- match.arg(metric)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  mu <- x$set_mu
  var <- x$set_var
  if (is.null(rownames(mu))) stop("set embedding must have rownames")

  all_sets <- rownames(mu)
  if (is.null(sets) && is.null(other_sets)) {
    return(gsemb_set_gaussian_distance(mu, var, metric = metric, eps = eps))
  }
  if (is.null(sets)) sets <- all_sets
  sets <- intersect(as.character(sets), all_sets)
  if (length(sets) == 0) stop("no sets found in embedding")
  if (is.null(other_sets)) {
    other_sets <- sets
  } else {
    other_sets <- intersect(as.character(other_sets), all_sets)
    if (length(other_sets) == 0) stop("no other_sets found in embedding")
  }
  gsemb_set_gaussian_distance(mu[sets, , drop = FALSE], var[sets, , drop = FALSE], mu[other_sets, , drop = FALSE], var[other_sets, , drop = FALSE], metric = metric, eps = eps)
}

#' Compute distances between two collections of gene clusters
#'
#' Gene clusters are represented as diagonal Gaussians fitted from member genes
#' in the learned gene embedding.
#'
#' @param x A `gsemb_embedding` object.
#' @param clusters Named list of character vectors (gene IDs).
#' @param other_clusters Optional second named list of clusters.
#' @param metric Distance metric: `"w2"` or `"sym_kl"`.
#' @param eps Lower bound for variances.
#' @param min_size Minimum number of genes per cluster after intersecting with the embedding.
#'
#' @return A numeric distance matrix (clusters x other_clusters).
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Define some custom clusters
#' clusters <- list(
#'   CLUST1 = c("GENE1", "GENE2", "GENE3"),
#'   CLUST2 = c("GENE4", "GENE5", "GENE6")
#' )
#' # Compute distances between these clusters
#' dist_mat <- gsemb_cluster_similarity(fit, clusters = clusters, metric = "w2")
#'
#' # Compute distances between two different cluster lists
#' other_clusters <- list(
#'   OTHER1 = c("GENE10", "GENE11"),
#'   OTHER2 = c("GENE12", "GENE13", "GENE14")
#' )
#' dist_mat2 <- gsemb_cluster_similarity(
#'   fit,
#'   clusters = clusters,
#'   other_clusters = other_clusters,
#'   metric = "sym_kl"
#' )
#' }
#' @export
gsemb_cluster_similarity <- function(x,
                                     clusters,
                                     other_clusters = NULL,
                                     metric = c("w2", "sym_kl"),
                                     eps = 1e-8,
                                     min_size = 2) {
  metric <- match.arg(metric)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  if (!is.list(clusters) || is.null(names(clusters))) stop("clusters must be a named list")
  if (is.null(other_clusters)) other_clusters <- clusters
  if (!is.list(other_clusters) || is.null(names(other_clusters))) stop("other_clusters must be a named list")

  gene_emb <- x$gene_embedding
  if (is.null(rownames(gene_emb))) stop("gene embedding must have rownames")

  to_gauss <- function(gs) {
    gs <- lapply(gs, function(v) intersect(as.character(v), rownames(gene_emb)))
    gs <- gs[vapply(gs, length, integer(1)) >= min_size]
    if (length(gs) == 0) stop("no clusters have enough genes in embedding")
    g <- gsemb_fit_set_gaussians_from_members(gene_emb, gs, eps = eps)
    list(mu = g$mu, var = g$var)
  }

  a <- to_gauss(clusters)
  b <- to_gauss(other_clusters)
  gsemb_set_gaussian_distance(a$mu, a$var, b$mu, b$var, metric = metric, eps = eps)
}

#' Create concise gene sets from a fitted embedding
#'
#' Convenience wrapper around `gsemb_make_concise_gene_sets` using the embedding
#' stored in a `gsemb_embedding` object.
#'
#' @param x A `gsemb_embedding` object.
#' @param gene_sets Named list of gene sets (used when restricting candidates).
#' @param ... Passed to `gsemb_make_concise_gene_sets`.
#'
#' @return A named list of concise gene sets.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Define original gene sets
#' orig_sets <- list(
#'   PATH1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
#'   PATH2 = c("GENE3", "GENE4", "GENE6", "GENE7")
#' )
#' # Create concise versions
#' concise <- gsemb_concise_gene_sets(fit, gene_sets = orig_sets, top_n = 3)
#'
#' # View concise sets
#' str(concise)
#' }
#' @export
gsemb_concise_gene_sets <- function(x,
                                    gene_sets,
                                    ...) {
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  gene_sets <- validate_gene_sets(gene_sets)
  gsemb_make_concise_gene_sets(
    gene_embedding = x$gene_embedding,
    set_mu = x$set_mu,
    set_var = x$set_var,
    gene_sets = gene_sets,
    ...
  )
}

#' One-click training from PPI edges and gene sets
#'
#' @param ppi_edges PPI edge list as a data.frame.
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param ... Passed to `gsemb_fit`.
#'
#' @return A `gsemb_embedding` object.
#' @examples
#' \dontrun{
#' # Simulate a small PPI edge list
#' set.seed(42)
#' nodes <- paste0("GENE", 1:50)
#' edges <- data.frame(
#'   node1 = sample(nodes, 200, replace = TRUE),
#'   node2 = sample(nodes, 200, replace = TRUE),
#'   weight = runif(200, 0.5, 1)
#' )
#' # Remove self-loops
#' edges <- edges[edges$node1 != edges$node2, ]
#'
#' # Simulate a few gene sets
#' gene_sets <- list(
#'   SET1 = sample(nodes, 10),
#'   SET2 = sample(nodes, 8),
#'   SET3 = sample(nodes, 12)
#' )
#'
#' # One-click training
#' fit <- gsemb_train_embedding_from_ppi_and_genesets(
#'   ppi_edges = edges,
#'   gene_sets = gene_sets,
#'   method = "svd",
#'   dim = 16,
#'   k = 20,
#'   alpha = 0.5
#' )
#'
#' # Inspect the result
#' print(fit)
#' }
#' @export
gsemb_train_embedding_from_ppi_and_genesets <- function(ppi_edges,
                                                        gene_sets,
                                                        ...) {
  gsemb_fit(ppi = ppi_edges, gene_sets = gene_sets, ...)
}

#' Compute all default similarity matrices from a fitted embedding
#'
#' @param x A `gsemb_embedding` object.
#' @param eps_gene Small constant for cosine similarity stability.
#' @param eps_set Lower bound for Gaussian variances.
#'
#' @return A list with `gene_gene`, `set_set_w2`, and `set_set_kl`.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Compute all similarity matrices
#' all_sims <- gsemb_calculate_all_similarities(fit)
#'
#' # Access individual matrices
#' gene_gene_sim <- all_sims$gene_gene
#' set_set_w2 <- all_sims$set_set_w2
#' set_set_kl <- all_sims$set_set_kl
#' }
#' @export
gsemb_calculate_all_similarities <- function(x,
                                             eps_gene = 1e-12,
                                             eps_set = 1e-8) {
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  list(
    gene_gene = gsemb_gene_similarity(x, eps = eps_gene),
    set_set_w2 = gsemb_set_similarity(x, metric = "w2", eps = eps_set),
    set_set_kl = gsemb_set_similarity(x, metric = "sym_kl", eps = eps_set)
  )
}

#' Create concise gene sets (high-level wrapper)
#'
#' @param x A `gsemb_embedding` object.
#' @param gene_sets Named list of gene sets.
#' @param ... Passed to `gsemb_concise_gene_sets`.
#'
#' @return A named list of concise gene sets.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Define original gene sets
#' orig_sets <- list(
#'   PATH1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
#'   PATH2 = c("GENE3", "GENE4", "GENE6", "GENE7")
#' )
#' # Get concise versions using the wrapper
#' concise <- gsemb_get_concise_gene_sets(fit, gene_sets = orig_sets, top_n = 3)
#'
#' # View concise sets
#' str(concise)
#' }
#' @export
gsemb_get_concise_gene_sets <- function(x,
                                        gene_sets,
                                        ...) {
  gsemb_concise_gene_sets(x, gene_sets = gene_sets, ...)
}

# 【更改前】置换零分布在 gsemb_embedding_enrichment 内 for 循环逐次 sample + crossprod（O(nperm) 次 R 调用）。
# 【更改后】优先 perm_mat + 一次 crossprod（BLAS）；其次 Rcpp；再次 vapply 生成 perm_mat 后 crossprod。
# perm_mat 由外层预生成，各批共用，保证分批时 p 值一致。
.gsemb_enrichment_null_scores <- function(W, stats_vec, nperm, seed, perm_mat = NULL, n_workers = 1L) {
  stats_vec <- as.numeric(stats_vec)

  if (!is.null(perm_mat)) {
    return(t(crossprod(W, perm_mat)))
  }

  if (.native_routine_available("_geneSetEmbedding_rcpp_enrichment_permutations")) {
    return(rcpp_enrichment_permutations(W, stats_vec, nperm, seed))
  }

  set.seed(seed)
  perm_inner <- vapply(
    seq_len(nperm),
    function(b) sample(stats_vec, length(stats_vec), replace = FALSE),
    FUN.VALUE = numeric(length(stats_vec))
  )

  if (n_workers <= 1L || !requireNamespace("future.apply", quietly = TRUE)) {
    return(t(crossprod(W, perm_inner)))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  null_list <- future.apply::future_lapply(
    seq_len(nperm),
    function(b) as.numeric(crossprod(W, perm_inner[, b])),
    future.seed = TRUE
  )
  do.call(rbind, null_list)
}

# 【更改前】gsemb_embedding_enrichment 为单函数体：一次算全部 sets 的 S/W/ES/置换/p.adjust。
# 【更改后】提取为单次 sets 批次的内部函数；ES 与置换共用 as.numeric(gene_stats) 顺序；
#          two.sided 用矩阵对矩阵比较；p.adjust 由外层合并后统一 BH。
.gsemb_embedding_enrichment_sets <- function(gene_stats,
                                             gene_emb,
                                             genes,
                                             mu,
                                             var,
                                             gene_sets,
                                             score,
                                             temperature,
                                             nperm,
                                             alternative,
                                             seed,
                                             eps,
                                             top_genes,
                                             perm_mat) {
  S <- gsemb_gene_to_set_score(
    gene_embedding = gene_emb,
    set_mu = mu,
    set_var = var,
    score = score,
    eps = eps
  )

  S <- as.matrix(S)
  col_max <- apply(S, 2, max)
  S <- S - rep(col_max, each = nrow(S))
  exp_S <- exp(S / temperature)
  col_sum <- colSums(exp_S)
  W <- exp_S / rep(col_sum, each = nrow(exp_S))
  W[!is.finite(W)] <- 0

  # 【更改前】crossprod(W, gene_stats) 使用 intersect() 顺序的命名向量，与置换 sample 顺序可能不一致。
  # 【更改后】as.numeric(gene_stats) 与 rownames(gene_emb) 严格同序，ES 与置换共用同一顺序。
  stats_vec <- as.numeric(gene_stats)
  es <- as.numeric(crossprod(W, stats_vec))
  names(es) <- colnames(W)
  n_es <- length(es)

  z <- rep(NA_real_, n_es)
  pvals <- rep(NA_real_, n_es)

  if (nperm > 0L) {
    null_scores <- .gsemb_enrichment_null_scores(W, stats_vec, nperm, seed, perm_mat)
    null_mean <- colMeans(null_scores)
    null_sd <- pmax(apply(null_scores, 2, stats::sd), eps)
    z <- (es - null_mean) / null_sd

    es_mat <- matrix(es, nrow = nperm, ncol = n_es, byrow = TRUE)
    if (alternative == "greater") {
      pvals <- colMeans(null_scores >= es_mat)
    } else if (alternative == "less") {
      pvals <- colMeans(null_scores <= es_mat)
    } else {
      # 【更改前】abs(t(t(null_scores) - cen)) >= abs(es - cen) 中矩阵与向量比较可能因回收算错。
      # 【更改后】cen_mat/shift_mat 均为 nperm x n_es，矩阵对矩阵比较。
      cen_mat <- matrix(null_mean, nrow = nperm, ncol = n_es, byrow = TRUE)
      shift_mat <- matrix(abs(es - null_mean), nrow = nperm, ncol = n_es, byrow = TRUE)
      pvals <- colMeans(abs(null_scores - cen_mat) >= shift_mat)
    }
    pvals <- pmax(pvals, 1 / nperm)
  }

  set_size <- rep(NA_integer_, n_es)
  if (!is.null(gene_sets)) {
    set_size <- vapply(names(es), function(sid) {
      if (!sid %in% names(gene_sets)) {
        return(NA_integer_)
      }
      length(intersect(gene_sets[[sid]], genes))
    }, integer(1))
  }

  core <- rep(NA_character_, length(es))
  if (top_genes > 0L) {
    for (j in seq_along(es)) {
      wj <- W[, j]
      ord <- order(wj, decreasing = TRUE)
      ord <- ord[seq_len(min(top_genes, length(ord)))]
      core[j] <- paste0(rownames(W)[ord], collapse = "/")
    }
  }

  data.frame(
    ID = names(es),
    ES = as.numeric(es),
    z = as.numeric(z),
    pvalue = as.numeric(pvals),
    p.adjust = NA_real_,
    setSize = as.integer(set_size),
    core_enrichment = core,
    stringsAsFactors = FALSE
  )
}

#' Enrichment analysis using gene‑set Gaussian embeddings
#'
#' Given a vector of gene‑level statistics (e.g., log‑fold‑changes, p‑values),
#' compute enrichment scores for each gene set based on the likelihood of the
#' statistics under the set’s Gaussian embedding. Permutation‑based p‑values
#' and multiple‑testing adjusted q‑values are provided.
#'
#' @param gene_stats Named numeric vector of gene‑level statistics.
#' @param x A `gsemb_embedding` object.
#' @param sets Optional character vector of set IDs to test.
#' @param gene_sets Named list of original gene‑set members (used for set‑size
#'   reporting and for restricting the candidate genes when `sets` is supplied).
#' @param score Scoring method: `"loglik"` (full Gaussian log‑likelihood) or
#'   `"neg_mahalanobis"` (negative squared Mahalanobis distance).
#' @param temperature Soft‑max temperature for converting scores to weights.
#' @param nperm Number of permutations for null distribution.
#' @param alternative Alternative hypothesis for permutation test.
#' @param seed Random seed for permutations.
#' @param eps Small constant added to variances for numerical stability.
#' @param top_genes Number of top‑weighted genes to report in `core_enrichment`.
#'
#' @return A data.frame with columns: `ID`, `ES` (enrichment score), `z` (z‑score),
#'   `pvalue`, `p.adjust`, `setSize`, `core_enrichment`.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Simulate some gene‑level statistics
#' all_genes <- rownames(fit$gene_embedding)
#' stats <- rnorm(length(all_genes))
#' names(stats) <- all_genes
#'
#' # Run enrichment analysis
#' enrich <- gsemb_embedding_enrichment(
#'   gene_stats = stats,
#'   x = fit,
#'   gene_sets = list(
#'     SET1 = sample(all_genes, 10),
#'     SET2 = sample(all_genes, 8)
#'   ),
#'   score = "loglik",
#'   nperm = 100,
#'   alternative = "greater"
#' )
#'
#' # View top enriched sets
#' head(enrich[order(enrich$pvalue), ])
#' }
#' @export
gsemb_embedding_enrichment <- function(gene_stats,
                                       x,
                                       sets = NULL,
                                       gene_sets = NULL,
                                       score = c("loglik", "neg_mahalanobis"),
                                       temperature = 1.0,
                                       nperm = 1000,
                                       alternative = c("two.sided", "greater", "less"),
                                       seed = 1,
                                       eps = 1e-8,
                                       top_genes = 30) {
  score <- match.arg(score)
  alternative <- match.arg(alternative)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  if (!is.numeric(gene_stats) || is.null(names(gene_stats))) stop("gene_stats must be a named numeric vector")
  if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) stop("temperature must be a positive scalar")
  if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 0) stop("nperm must be a non-negative integer")
  nperm <- as.integer(nperm)
  if (!is.numeric(top_genes) || length(top_genes) != 1 || top_genes < 0) stop("top_genes must be a non-negative integer")
  top_genes <- as.integer(top_genes)

  gene_emb <- x$gene_embedding
  if (is.null(rownames(gene_emb))) stop("gene embedding must have rownames")
  mu <- as_numeric_matrix(x$set_mu)
  var <- as_numeric_matrix(x$set_var)
  if (is.null(rownames(mu)) || is.null(rownames(var))) stop("set embedding must have rownames")
  if (nrow(mu) != nrow(var)) stop("set_mu and set_var must have the same number of rows")
  if (length(rownames(mu)) != nrow(mu)) {
    stop("set_mu rownames length (", length(rownames(mu)), ") != nrow (", nrow(mu), ")")
  }

  genes <- intersect(names(gene_stats), rownames(gene_emb))
  if (length(genes) < 2) stop("not enough genes overlap between gene_stats and gene_embedding")
  gene_emb <- gene_emb[genes, , drop = FALSE]
  # 【更改前】gene_stats <- gene_stats[genes]，顺序为 intersect() 返回序（非 embedding 行序）。
  # 【更改后】按 rownames(gene_emb) 重排，与 ES/置换的 as.numeric 顺序严格一致。
  gene_stats <- gene_stats[rownames(gene_emb)]
  storage.mode(gene_stats) <- "double"

  set_names <- rownames(mu)
  if (is.null(sets)) {
    set_idx <- seq_len(nrow(mu))
  } else {
    set_idx <- match(as.character(sets), set_names, nomatch = 0L)
    if (any(set_idx == 0L)) {
      bad <- sets[set_idx == 0L]
      stop(
        "some sets are not in set_mu rownames; first missing: ",
        paste(head(bad, 3L), collapse = ", ")
      )
    }
  }
  if (length(set_idx) == 0L) stop("no sets found in embedding")
  mu <- mu[set_idx, , drop = FALSE]
  var <- var[set_idx, , drop = FALSE]
  sets <- set_names[set_idx]

  if (!is.null(gene_sets)) {
    gene_sets <- validate_gene_sets(gene_sets)
  }

  # 【更改前】无预生成 perm_mat；置换在函数内 set.seed 后逐次 sample（分批时会各自生成）。
  # 【更改后】外层 set.seed 一次生成 perm_mat，各批 .gsemb_embedding_enrichment_sets 共用。
  perm_mat <- NULL
  if (nperm > 0L) {
    stats_vec <- as.numeric(gene_stats)
    set.seed(seed)
    perm_mat <- vapply(
      seq_len(nperm),
      function(b) sample(stats_vec, length(stats_vec), replace = FALSE),
      FUN.VALUE = numeric(length(stats_vec))
    )
  }

  # 【更改前】一次 gsemb_gene_to_set_score(全部 sets)，n_g*n_s 过大时 OOM。
  # 【更改后】n_g*n_s > 2e7 时按 sets 自动拆批，rbind 后统一 p.adjust(BH)。
  n_sets <- length(set_idx)
  n_genes <- nrow(gene_emb)
  sets_per_batch <- if (as.numeric(n_genes) * as.numeric(n_sets) > 2e7) {
    max(50L, as.integer(floor(2e7 / max(n_genes, 1L))))
  } else {
    n_sets
  }
  if (!is.finite(sets_per_batch) || sets_per_batch < 1L) {
    sets_per_batch <- n_sets
  }

  idx_batches <- split(set_idx, ceiling(seq_along(set_idx) / sets_per_batch))
  parts <- lapply(idx_batches, function(batch_idx) {
    if (any(batch_idx < 1L | batch_idx > nrow(mu))) {
      stop("invalid set index in batch: min=", min(batch_idx), ", max=", max(batch_idx), ", nrow=", nrow(mu))
    }
    .gsemb_embedding_enrichment_sets(
      gene_stats = gene_stats,
      gene_emb = gene_emb,
      genes = genes,
      mu = mu[batch_idx, , drop = FALSE],
      var = var[batch_idx, , drop = FALSE],
      gene_sets = gene_sets,
      score = score,
      temperature = temperature,
      nperm = nperm,
      alternative = alternative,
      seed = seed,
      eps = eps,
      top_genes = top_genes,
      perm_mat = perm_mat
    )
  })

  out <- do.call(rbind, parts)
  out$p.adjust <- stats::p.adjust(out$pvalue, method = "BH")
  out
}
