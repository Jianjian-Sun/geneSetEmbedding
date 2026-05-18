.rcpp_available <- function() {
  requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)
}

#' Compute cosine similarity between gene embeddings
#'
#' @param gene_embedding Numeric matrix with genes in rows and embedding dimensions in columns.
#'   Row names must be gene identifiers.
#' @param other Optional numeric matrix with the same number of columns as `gene_embedding`.
#'   If `NULL`, similarities are computed within `gene_embedding`.
#' @param eps Numeric scalar added to norms for numerical stability.
#'
#' @return A numeric matrix of cosine similarities. Row/column names are inherited from
#'   the row names of the input embedding matrices.
#' @examples
#' # Create two random embedding matrices
#' set.seed(42)
#' X <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' rownames(X) <- paste0("GENE", 1:5)
#' Y <- matrix(rnorm(12), nrow = 3, ncol = 4)
#' rownames(Y) <- paste0("OTHER", 1:3)
#'
#' # Compute cosine similarity within X
#' sim1 <- gsemb_gene_cosine_similarity(X)
#' dim(sim1)
#'
#' # Compute cosine similarity between X and Y
#' sim2 <- gsemb_gene_cosine_similarity(X, Y)
#' dim(sim2)
#' @export
gsemb_gene_cosine_similarity <- function(gene_embedding, other = NULL, eps = 1e-12) {
  X <- as_numeric_matrix(gene_embedding)
  if (is.null(rownames(X))) stop("gene_embedding must have rownames")
  if (is.null(other)) {
    Y <- X
  } else {
    Y <- as_numeric_matrix(other)
    if (is.null(rownames(Y))) stop("other must have rownames")
  }
  if (ncol(X) != ncol(Y)) stop("embeddings must have the same number of columns")
  x_norm <- sqrt(rowSums(X * X))
  y_norm <- sqrt(rowSums(Y * Y))
  x_norm <- pmax(x_norm, eps)
  y_norm <- pmax(y_norm, eps)
  sim <- (X / x_norm) %*% t(Y / y_norm)
  rownames(sim) <- rownames(X)
  colnames(sim) <- rownames(Y)
  sim
}

#' Compute pairwise distances between diagonal Gaussian set embeddings
#'
#' @param set_mu Numeric matrix of Gaussian means (gene sets in rows).
#' @param set_var Numeric matrix of diagonal variances (same shape/rownames as `set_mu`).
#' @param other_mu Optional matrix of means for the second collection.
#' @param other_var Optional matrix of variances for the second collection.
#' @param metric Distance metric: `"w2"` for diagonal 2-Wasserstein distance, or `"sym_kl"`
#'   for symmetric KL divergence.
#' @param eps Numeric scalar lower bound for variances.
#'
#' @return A numeric matrix of pairwise distances.
#' @examples
#' # Create example Gaussian parameters
#' set.seed(123)
#' mu1 <- matrix(rnorm(6), nrow = 2, ncol = 3)
#' var1 <- matrix(rexp(6, rate = 2), nrow = 2, ncol = 3)
#' rownames(mu1) <- rownames(var1) <- c("SET_A", "SET_B")
#'
#' mu2 <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' var2 <- matrix(rexp(9, rate = 2), nrow = 3, ncol = 3)
#' rownames(mu2) <- rownames(var2) <- c("SET_X", "SET_Y", "SET_Z")
#'
#' # Compute Wasserstein-2 distances between the two collections
#' dist_w2 <- gsemb_set_gaussian_distance(mu1, var1, mu2, var2, metric = "w2")
#' dim(dist_w2)
#'
#' # Compute symmetric KL divergences
#' dist_kl <- gsemb_set_gaussian_distance(mu1, var1, mu2, var2, metric = "sym_kl")
#' dim(dist_kl)
#' @export
gsemb_set_gaussian_distance <- function(set_mu,
                                        set_var,
                                        other_mu = NULL,
                                        other_var = NULL,
                                        metric = c("w2", "sym_kl"),
                                        eps = 1e-8) {
  metric <- match.arg(metric)
  mu <- as_numeric_matrix(set_mu)
  var <- as_numeric_matrix(set_var)
  if (is.null(rownames(mu)) || is.null(rownames(var))) stop("set_mu/set_var must have rownames")
  if (!all(dim(mu) == dim(var))) stop("set_mu and set_var must have identical dimensions")
  if (!all(rownames(mu) == rownames(var))) stop("set_mu and set_var must have identical rownames")

  if (is.null(other_mu)) other_mu <- mu
  if (is.null(other_var)) other_var <- var
  mu2 <- as_numeric_matrix(other_mu)
  var2 <- as_numeric_matrix(other_var)
  if (!all(dim(mu2) == dim(var2))) stop("other_mu and other_var must have identical dimensions")
  if (ncol(mu) != ncol(mu2)) stop("mu and other_mu must have same number of columns")

  var <- pmax(var, eps)
  var2 <- pmax(var2, eps)

  if (.rcpp_available()) {
    cpp_out <- if (metric == "w2") {
      w2_distance(mu, var, mu2, var2)
    } else {
      sym_kl_distance(mu, var, mu2, var2)
    }
    rownames(cpp_out) <- rownames(mu)
    colnames(cpp_out) <- rownames(mu2)
    return(cpp_out)
  }

  # Fallback: pure-R vectorized implementation
  if (metric == "w2") {
    mu_sq <- rowSums(mu * mu)
    mu2_sq <- rowSums(mu2 * mu2)
    dmu2 <- mu_sq + outer(mu2_sq, mu_sq, "+") - 2 * mu %*% t(mu2)
    sv <- sqrt(var)
    sv2 <- sqrt(var2)
    var_term <- rowSums(sv * sv) + outer(rowSums(sv2 * sv2), rep(1, nrow(mu)), "+") - 2 * sv %*% t(sv2)
    out <- dmu2 + var_term
  } else {
    inv_var <- 1 / var
    inv_var2 <- 1 / var2
    ratio_mat <- var %*% t(inv_var2)
    ratio_mat2 <- var2 %*% t(inv_var)
    mahal_ij <- (
      t(t(rowSums(mu * mu)) %*% t(inv_var2)) +
        t(t(rowSums(mu2 * mu2)) %*% t(inv_var2)) -
        2 * mu %*% (t(inv_var2) %*% mu2)
    )
    mahal_ji <- (
      t(t(rowSums(mu2 * mu2)) %*% t(inv_var)) +
        t(t(rowSums(mu * mu)) %*% t(inv_var)) -
        2 * mu2 %*% (t(inv_var) %*% mu)
    )
    out <- 0.5 * (ratio_mat + t(ratio_mat2) + mahal_ij + mahal_ji)
  }

  rownames(out) <- rownames(mu)
  colnames(out) <- rownames(mu2)
  out
}
