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

  inv_var <- 1 / var # n_s x d

  # Tile gene matrix (n_g x d) so each row of E is repeated n_s times:
  # E_tiled: (n_g * n_s) x d  [row order: set1_gene1, set1_gene2, ..., set2_gene1, ...]
  E_tiled <- E[rep(seq_len(nrow(E)), each = nrow(mu)), , drop = FALSE]
  # Tile mu (n_s x d) so each set mean is repeated n_g times, same block order:
  # mu_tiled: (n_g * n_s) x d
  mu_tiled <- mu[rep(seq_len(nrow(mu)), times = nrow(E)), , drop = FALSE]
  # Tile inv_var (n_s x d) the same way:
  inv_var_tiled <- inv_var[rep(seq_len(nrow(mu)), times = nrow(E)), , drop = FALSE]

  # Squared differences weighted by inverse variance: (n_g * n_s) x d
  diff_sq <- (E_tiled - mu_tiled)^2 * inv_var_tiled
  # Sum over dimensions -> weighted Mahalanobis distances: (n_g * n_s)
  md2 <- rowSums(diff_sq)

  # Reshape back to n_g x n_s matrix
  out <- matrix(md2, nrow = nrow(E), ncol = nrow(mu), byrow = TRUE)

  if (score == "neg_mahalanobis") {
    out <- -out
  } else {
    logdet <- rowSums(log(var)) # n_s
    out <- -0.5 * (out + rep(logdet, each = nrow(E)))
  }
  rownames(out) <- rownames(E)
  colnames(out) <- rownames(mu)
  out
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
      x <- (s_sorted - max(s_sorted)) / temperature
      w <- exp(x)
      w <- w / sum(w)
      cumw <- cumsum(w)
      k <- which(cumw >= mass)[1]
      if (is.na(k) || k < 1) k <- 1
      if (k < min_size) k <- min_size
      if (k > max_size) k <- max_size
      k <- min(k, length(candidates))
      concise[[sid]] <- candidates[seq_len(k)]
    }
  }
  concise <- concise[vapply(concise, length, integer(1)) > 0]
  concise
}
