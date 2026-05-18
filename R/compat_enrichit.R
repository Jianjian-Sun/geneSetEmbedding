#' Coerce an embedding object to a numeric matrix
#'
#' @param x A numeric matrix, Matrix, or data.frame representing an embedding.
#' @param id_col Optional ID column name when `x` is a data.frame.
#'
#' @return A numeric matrix.
#' @examples
#' # Create a data.frame with an ID column
#' df <- data.frame(
#'   gene = paste0("GENE", 1:10),
#'   dim1 = rnorm(10),
#'   dim2 = rnorm(10)
#' )
#' # Convert to numeric matrix using the 'gene' column as row names
#' mat <- gsemb_as_embedding_matrix(df, id_col = "gene")
#' dim(mat)
#' rownames(mat)
#' @export
gsemb_as_embedding_matrix <- function(x, id_col = NULL) {
  as_numeric_matrix(x, id_col = id_col)
}

#' Row-wise cosine similarity
#'
#' Compatibility wrapper for cosine similarity.
#'
#' @param x Numeric embedding matrix.
#' @param y Optional second embedding matrix.
#' @param eps Small constant for numerical stability.
#'
#' @return A cosine similarity matrix.
#' @examples
#' # Create two random embedding matrices
#' set.seed(42)
#' X <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' rownames(X) <- paste0("ROW", 1:5)
#' Y <- matrix(rnorm(12), nrow = 3, ncol = 4)
#' rownames(Y) <- paste0("OTHER", 1:3)
#'
#' # Compute cosine similarity within X
#' sim1 <- gsemb_row_cosine_similarity(X)
#' dim(sim1)
#'
#' # Compute cosine similarity between X and Y
#' sim2 <- gsemb_row_cosine_similarity(X, Y)
#' dim(sim2)
#' @export
gsemb_row_cosine_similarity <- function(x, y = NULL, eps = 1e-12) {
  gsemb_gene_cosine_similarity(gene_embedding = x, other = y, eps = eps)
}

#' Diagonal Gaussian 2-Wasserstein distance
#'
#' Compatibility wrapper for diagonal Gaussian distance with `metric="w2"`.
#'
#' @param mu Mean matrix for sets.
#' @param var Variance matrix for sets.
#' @param mu2 Optional mean matrix for other sets.
#' @param var2 Optional variance matrix for other sets.
#' @param eps Lower bound for variances.
#'
#' @return A distance matrix.
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
#' dist_w2 <- gsemb_diag_gaussian_w2(mu1, var1, mu2, var2)
#' dim(dist_w2)
#' @export
gsemb_diag_gaussian_w2 <- function(mu, var, mu2 = NULL, var2 = NULL, eps = 1e-8) {
  gsemb_set_gaussian_distance(set_mu = mu, set_var = var, other_mu = mu2, other_var = var2, metric = "w2", eps = eps)
}

#' Diagonal Gaussian symmetric KL distance
#'
#' Compatibility wrapper for diagonal Gaussian distance with `metric="sym_kl"`.
#'
#' @param mu Mean matrix for sets.
#' @param var Variance matrix for sets.
#' @param mu2 Optional mean matrix for other sets.
#' @param var2 Optional variance matrix for other sets.
#' @param eps Lower bound for variances.
#'
#' @return A distance matrix.
#' @examples
#' # Use the same example data as in gsemb_diag_gaussian_w2
#' set.seed(123)
#' mu1 <- matrix(rnorm(6), nrow = 2, ncol = 3)
#' var1 <- matrix(rexp(6, rate = 2), nrow = 2, ncol = 3)
#' rownames(mu1) <- rownames(var1) <- c("SET_A", "SET_B")
#'
#' mu2 <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' var2 <- matrix(rexp(9, rate = 2), nrow = 3, ncol = 3)
#' rownames(mu2) <- rownames(var2) <- c("SET_X", "SET_Y", "SET_Z")
#'
#' # Compute symmetric KL divergences
#' dist_kl <- gsemb_diag_gaussian_sym_kl(mu1, var1, mu2, var2)
#' dim(dist_kl)
#' @export
gsemb_diag_gaussian_sym_kl <- function(mu, var, mu2 = NULL, var2 = NULL, eps = 1e-8) {
  gsemb_set_gaussian_distance(set_mu = mu, set_var = var, other_mu = mu2, other_var = var2, metric = "sym_kl", eps = eps)
}

#' Gene-to-set score under diagonal Gaussian
#'
#' Compatibility wrapper for `gsemb_gene_to_set_score`.
#'
#' @param gene_emb Numeric gene embedding matrix.
#' @param set_mu Mean matrix for sets.
#' @param set_var Variance matrix for sets.
#' @param score Scoring function.
#' @param eps Lower bound for variances.
#'
#' @return A numeric matrix of scores.
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
#' scores_loglik <- gsemb_gene_to_diag_gaussian_score(
#'   gene_emb, set_mu, set_var,
#'   score = "loglik"
#' )
#' dim(scores_loglik)
#'
#' # Compute negative Mahalanobis distance scores
#' scores_negmah <- gsemb_gene_to_diag_gaussian_score(
#'   gene_emb, set_mu, set_var,
#'   score = "neg_mahalanobis"
#' )
#' dim(scores_negmah)
#' @export
gsemb_gene_to_diag_gaussian_score <- function(gene_emb,
                                              set_mu,
                                              set_var,
                                              score = c("loglik", "neg_mahalanobis"),
                                              eps = 1e-8) {
  gsemb_gene_to_set_score(gene_embedding = gene_emb, set_mu = set_mu, set_var = set_var, score = score, eps = eps)
}

#' Make concise gene sets from diagonal Gaussian embeddings
#'
#' Compatibility wrapper for `gsemb_make_concise_gene_sets`.
#'
#' @param gene_emb Numeric gene embedding matrix.
#' @param set_mu Mean matrix for sets.
#' @param set_var Variance matrix for sets.
#' @param gene_sets Optional named list of gene sets used to restrict candidates.
#' @param top_n Number of genes per set when using `select="top_n"`.
#' @param restrict_to_members Restrict candidate genes to members of the input set.
#' @param min_size Minimum size of each output set.
#' @param max_size Maximum size of each output set.
#' @param score Scoring function.
#' @param eps Lower bound for variances.
#' @param ... Passed to `gsemb_make_concise_gene_sets`.
#'
#' @return A named list of concise gene sets.
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
#' concise <- gsemb_make_concise_gene_sets_from_gaussian(
#'   gene_emb, set_mu, set_var,
#'   gene_sets = gene_sets,
#'   top_n = 3,
#'   restrict_to_members = TRUE,
#'   score = "loglik"
#' )
#' str(concise)
#' @export
gsemb_make_concise_gene_sets_from_gaussian <- function(gene_emb,
                                                       set_mu,
                                                       set_var,
                                                       gene_sets = NULL,
                                                       top_n = 50,
                                                       restrict_to_members = TRUE,
                                                       min_size = 5,
                                                       max_size = 500,
                                                       score = c("loglik", "neg_mahalanobis"),
                                                       eps = 1e-8,
                                                       ...) {
  gsemb_make_concise_gene_sets(
    gene_embedding = gene_emb,
    set_mu = set_mu,
    set_var = set_var,
    gene_sets = gene_sets,
    top_n = top_n,
    restrict_to_members = restrict_to_members,
    min_size = min_size,
    max_size = max_size,
    score = score,
    eps = eps,
    ...
  )
}
