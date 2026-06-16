#' Score genes against diagonal Gaussian set embeddings
#'
#' Compute gene-to-set scores under a diagonal Gaussian representation of sets.
#' This is used for selecting concise gene sets.
#'
#' @param gene_embedding Numeric matrix of gene embeddings (genes in rows).
#' @param set_mu Numeric matrix of set means (sets in rows).
#' @param set_var Numeric matrix of set diagonal variances (same shape as `set_mu`).
#' @param score Scoring function: `"loglik"` or `"neg_mahalanobis"`.
#' @param eps Lower bound for variances.
#'
#' @return A numeric matrix with genes in rows and sets in columns.
#' @examples
#' # Simulate a small gene embedding and set parameters
#' set.seed(456)
#' gene_emb <- matrix(rnorm(15), nrow = 5, ncol = 3)
#' rownames(gene_emb) <- paste0("GENE", 1:5)
#'
#' set_mu <- matrix(rnorm(6), nrow = 2, ncol = 3)
#' set_var <- matrix(rexp(6, rate = 1), nrow = 2, ncol = 3)
#' rownames(set_mu) <- rownames(set_var) <- c("PATH1", "PATH2")
#'
#' # Compute log-likelihood scores
#' scores_loglik <- gsemb_gene_to_set_score(
#'   gene_emb, set_mu, set_var,
#'   score = "loglik"
#' )
#' dim(scores_loglik)
#'
#' # Compute negative Mahalanobis distance scores
#' scores_negmah <- gsemb_gene_to_set_score(
#'   gene_emb, set_mu, set_var,
#'   score = "neg_mahalanobis"
#' )
#' dim(scores_negmah)
#' @export
gsemb_gene_to_set_score <- function(gene_embedding,
                                    set_mu,
                                    set_var,
                                    score = c("loglik", "neg_mahalanobis"),
                                    eps = 1e-8) {
  score <- match.arg(score)
  E <- as_numeric_matrix(gene_embedding)
  mu <- as_numeric_matrix(set_mu)
  var <- as_numeric_matrix(set_var)
  if (is.null(rownames(E)) || is.null(rownames(mu))) stop("embeddings must have rownames")
  if (!all(dim(mu) == dim(var))) stop("set_mu and set_var must have the same dimensions")
  if (ncol(E) != ncol(mu)) stop("gene_embedding and set_mu must have the same number of columns")
  var <- pmax(var, eps)
  inv_var <- 1 / var

  # 【更改前】将 E/mu/inv_var 扩成 (n_g*n_s) x d 的巨型临时矩阵，再 rowSums 后 reshape；
  #          峰值内存约 n_g*n_s*d（Reactome 规模可达 8–9 GB）。
  # 【更改后】用 BLAS 矩阵乘法直接得到 n_g x n_s 的 md2 矩阵，峰值约 n_g*n_s（几百 MB）：
  #   md2[g,s] = sum_d (E[g,d]-mu[s,d])^2/var[s,d]
  #            = E^2 %*% t(inv_var) - 2*E %*% t(mu*inv_var) + sum_d mu^2/var
  E2 <- E * E
  mu_w <- mu * inv_var
  md2 <- E2 %*% t(inv_var) - 2 * (E %*% t(mu_w))
  md2 <- sweep(md2, 2L, rowSums(mu * mu * inv_var), "+")

  if (score == "neg_mahalanobis") {
    out <- -md2
  } else {
    logdet <- rowSums(log(var))
    out <- -0.5 * sweep(md2, 2L, logdet, "+")
  }
  rownames(out) <- rownames(E)
  colnames(out) <- rownames(mu)
  out
}

# 【更改前】gsemb_make_concise_gene_sets 在 select="softmax_mass" 时对每个基因集单独 exp + cumsum。
# 【更改后】restrict_to_members=FALSE 时用 .gsemb_softmax_mass_k_matrix 对所有列一次向量化；
#          单集路径仍用 .gsemb_softmax_mass_k，逻辑与原先循环体相同。
.gsemb_softmax_mass_k <- function(s_sorted, mass, temperature, min_size, max_size) {
  x <- (s_sorted - max(s_sorted)) / temperature
  w <- exp(x)
  w <- w / sum(w)
  cumw <- cumsum(w)
  k <- which(cumw >= mass)[1]
  if (is.na(k) || k < 1) k <- 1
  k <- max(k, min_size)
  k <- min(k, max_size, length(s_sorted))
  k
}

.gsemb_softmax_mass_k_matrix <- function(scores_mat, mass, temperature, min_size, max_size) {
  S <- scores_mat - rep(apply(scores_mat, 2, max), each = nrow(scores_mat))
  W <- exp(S / temperature)
  W <- W / rep(colSums(W), each = nrow(W))
  cumw <- apply(W, 2, cumsum)
  k_raw <- apply(cumw >= mass, 2, function(col) which(col)[1])
  k_raw[is.na(k_raw)] <- 1L
  pmax(min_size, pmin(max_size, k_raw, nrow(scores_mat)))
}

#' Create concise gene sets from embeddings
#'
#' Rank candidate genes by their score to each set's Gaussian embedding and
#' select a subset according to a selection rule.
#'
#' @param gene_embedding Numeric matrix of gene embeddings (genes in rows).
#' @param set_mu Numeric matrix of set means (sets in rows).
#' @param set_var Numeric matrix of set diagonal variances (same shape as `set_mu`).
#' @param gene_sets Optional named list of character vectors, used when restricting
#'   candidates to original members.
#' @param top_n Number of genes when `select="top_n"`.
#' @param select Selection strategy: `"top_n"`, `"softmax_mass"`, `"quantile"`, or `"score_threshold"`.
#' @param mass Softmax cumulative mass to cover when `select="softmax_mass"`.
#' @param temperature Softmax temperature when `select="softmax_mass"`.
#' @param quantile Quantile threshold when `select="quantile"`.
#' @param score_threshold Score cutoff when `select="score_threshold"`.
#' @param restrict_to_members Logical; restrict candidate genes to `gene_sets[[set]]`.
#' @param min_size Minimum size of each output set.
#' @param max_size Maximum size of each output set.
#' @param score Scoring function passed to `gsemb_gene_to_set_score`.
#' @param eps Lower bound for variances.
#'
#' @return A named list of concise gene sets (character vectors).
#' @examples
#' # Simulate embeddings and gene sets
#' set.seed(789)
#' gene_emb <- matrix(rnorm(30), nrow = 10, ncol = 3)
#' rownames(gene_emb) <- paste0("GENE", 1:10)
#'
#' set_mu <- matrix(rnorm(6), nrow = 2, ncol = 3)
#' set_var <- matrix(rexp(6, rate = 0.5), nrow = 2, ncol = 3)
#' rownames(set_mu) <- rownames(set_var) <- c("PATHWAY_A", "PATHWAY_B")
#'
#' gene_sets <- list(
#'   PATHWAY_A = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
#'   PATHWAY_B = c("GENE6", "GENE7", "GENE8", "GENE9", "GENE10")
#' )
#'
#' # Create concise gene sets (top 3 genes per set)
#' concise <- gsemb_make_concise_gene_sets(
#'   gene_emb, set_mu, set_var,
#'   gene_sets = gene_sets,
#'   top_n = 3,
#'   restrict_to_members = TRUE,
#'   score = "loglik"
#' )
#' str(concise)
#' @export
gsemb_make_concise_gene_sets <- function(gene_embedding,
                                         set_mu,
                                         set_var,
                                         gene_sets = NULL,
                                         top_n = 50,
                                         select = c("top_n", "softmax_mass", "quantile", "score_threshold"),
                                         mass = 0.8,
                                         temperature = 1.0,
                                         quantile = 0.9,
                                         score_threshold = NULL,
                                         restrict_to_members = TRUE,
                                         min_size = 5,
                                         max_size = 500,
                                         score = c("loglik", "neg_mahalanobis"),
                                         eps = 1e-8) {
  select <- match.arg(select)
  score <- match.arg(score)
  E <- as_numeric_matrix(gene_embedding)
  mu <- as_numeric_matrix(set_mu)
  var <- as_numeric_matrix(set_var)
  if (is.null(rownames(E))) stop("gene_embedding must have rownames")
  if (is.null(rownames(mu)) || is.null(rownames(var))) stop("set_mu/set_var must have rownames")
  if (!all(rownames(mu) == rownames(var))) stop("set_mu and set_var must have identical rownames")

  if (!is.null(gene_sets)) {
    gene_sets <- validate_gene_sets(gene_sets)
  } else {
    restrict_to_members <- FALSE
  }

  set_ids <- rownames(mu)
  concise <- vector("list", length(set_ids))
  names(concise) <- set_ids

  if (restrict_to_members && !is.null(gene_sets)) {
    # All sets share the same candidate pool per set -> compute all at once
    # to leverage the vectorized matrix product.
    all_scores <- gsemb_gene_to_set_score(E, mu, var, score = score, eps = eps)
  }

  # 【更改前】所有 select 模式均走下方 for 循环；softmax_mass 每集单独 exp/cumsum。
  # 【更改后】restrict_to_members=FALSE 且 softmax_mass 时，先算全矩阵再 .gsemb_softmax_mass_k_matrix 向量化。
  if (select == "softmax_mass" && !restrict_to_members) {
    if (!is.numeric(mass) || length(mass) != 1 || mass <= 0 || mass > 1) stop("mass must be in (0, 1]")
    if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) stop("temperature must be > 0")
    scores_mat <- gsemb_gene_to_set_score(E, mu, var, score = score, eps = eps)
    ord <- apply(scores_mat, 2, order, decreasing = TRUE)
    sorted_scores <- sapply(seq_len(ncol(scores_mat)), function(j) scores_mat[ord[, j], j])
    k_vec <- .gsemb_softmax_mass_k_matrix(sorted_scores, mass, temperature, min_size, max_size)
    for (j in seq_along(set_ids)) {
      sid <- set_ids[j]
      concise[[sid]] <- rownames(E)[ord[, j]][seq_len(k_vec[j])]
    }
  } else {
  for (sid in set_ids) {
    if (restrict_to_members) {
      if (!sid %in% names(gene_sets)) next
      candidates <- intersect(gene_sets[[sid]], rownames(E))
    } else {
      candidates <- rownames(E)
    }
    if (length(candidates) == 0) next
    s <- if (restrict_to_members && !is.null(gene_sets)) {
      all_scores[candidates, sid]
    } else {
      gsemb_gene_to_set_score(E[candidates, , drop = FALSE], mu[sid, , drop = FALSE], var[sid, , drop = FALSE], score = score, eps = eps)[, 1]
    }
    ord <- order(s, decreasing = TRUE)
    candidates <- candidates[ord]
    s_sorted <- s[ord]

    if (select == "top_n") {
      k <- min(top_n, length(candidates), max_size)
      k <- max(k, min_size)
      k <- min(k, length(candidates))
      concise[[sid]] <- candidates[seq_len(k)]
    } else if (select == "quantile") {
      if (is.na(quantile) || quantile <= 0 || quantile >= 1) stop("quantile must be in (0, 1)")
      thr <- as.numeric(stats::quantile(s_sorted, probs = quantile, type = 7, names = FALSE))
      keep <- which(s_sorted >= thr)
      k <- length(keep)
      if (k < min_size) k <- min_size
      if (k > max_size) k <- max_size
      k <- min(k, length(candidates))
      concise[[sid]] <- candidates[seq_len(k)]
    } else if (select == "score_threshold") {
      if (is.null(score_threshold) || !is.numeric(score_threshold) || length(score_threshold) != 1) {
        stop("score_threshold must be a numeric scalar when select='score_threshold'")
      }
      keep <- which(s_sorted >= score_threshold)
      k <- length(keep)
      if (k < min_size) k <- min_size
      if (k > max_size) k <- max_size
      k <- min(k, length(candidates))
      concise[[sid]] <- candidates[seq_len(k)]
    } else {
      if (!is.numeric(mass) || length(mass) != 1 || mass <= 0 || mass > 1) stop("mass must be in (0, 1]")
      if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) stop("temperature must be > 0")
      # 【更改前】内联 exp + cumsum 循环体（见 git 历史 concise.R 第 213–222 行）。
      # 【更改后】提取为 .gsemb_softmax_mass_k，语义不变。
      k <- .gsemb_softmax_mass_k(s_sorted, mass, temperature, min_size, max_size)
      concise[[sid]] <- candidates[seq_len(k)]
    }
  }
  }
  concise <- concise[vapply(concise, length, integer(1)) > 0]
  concise
}
