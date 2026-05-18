#' Build a PPI graph adjacency matrix
#'
#' Construct a sparse adjacency matrix from an edge list. Node names are stored
#' in the matrix dimnames and are used throughout the package.
#'
#' @param edges A data.frame containing at least two columns for endpoints.
#' @param node1,node2 Column names in \code{edges} for source/target.
#' @param weight Optional column name in \code{edges} providing edge weights.
#' @param directed Logical; whether to keep the graph directed.
#' @param nodes Optional character vector of node IDs to keep/order.
#'
#' @return A sparse adjacency matrix of class \code{dgCMatrix}.
#' @examples
#' # Create a small edge list
#' edges <- data.frame(
#'   node1 = c("A", "B", "C", "A"),
#'   node2 = c("B", "C", "A", "C"),
#'   weight = c(1.0, 2.0, 0.5, 1.5)
#' )
#'
#' # Build undirected weighted graph
#' adj <- gsemb_build_graph(edges, weight = "weight")
#' adj
#'
#' # Build directed unweighted graph
#' adj_dir <- gsemb_build_graph(edges, directed = TRUE)
#' adj_dir
#' @export
gsemb_build_graph <- function(edges,
                              node1 = "node1",
                              node2 = "node2",
                              weight = NULL,
                              directed = FALSE,
                              nodes = NULL) {
  if (!is.data.frame(edges)) stop("edges must be a data.frame")
  if (!node1 %in% names(edges)) stop("node1 column not found")
  if (!node2 %in% names(edges)) stop("node2 column not found")

  src <- as.character(edges[[node1]])
  dst <- as.character(edges[[node2]])
  if (is.null(weight)) {
    w <- rep(1, length(src))
  } else {
    if (!weight %in% names(edges)) stop("weight column not found")
    w <- as.numeric(edges[[weight]])
    w[is.na(w)] <- 0
  }

  if (is.null(nodes)) {
    nodes <- sort(unique(c(src, dst)))
  } else {
    nodes <- unique(as.character(nodes))
  }
  idx1 <- match(src, nodes)
  idx2 <- match(dst, nodes)
  ok <- !is.na(idx1) & !is.na(idx2) & w != 0
  idx1 <- idx1[ok]
  idx2 <- idx2[ok]
  w <- w[ok]

  adj <- Matrix::sparseMatrix(i = idx1, j = idx2, x = w, dims = c(length(nodes), length(nodes)))
  if (!directed) {
    adj <- adj + Matrix::t(adj)
  }
  dimnames(adj) <- list(nodes, nodes)
  Matrix::drop0(adj)
}

#' Construct a random-walk transition matrix
#'
#' Normalize an adjacency matrix into a transition matrix by column- or
#' row-normalization.
#'
#' @param adj A sparse/dense matrix (typically produced by \code{gsemb_build_graph}).
#' @param normalize One of \code{"col"} (column-stochastic) or \code{"row"} (row-stochastic).
#' @param eps Small constant to avoid division by zero.
#'
#' @return A sparse transition matrix.
#' @examples
#' # Build a small adjacency matrix
#' edges <- data.frame(
#'   node1 = c("A", "B", "C"),
#'   node2 = c("B", "C", "A")
#' )
#' adj <- gsemb_build_graph(edges)
#'
#' # Column-normalized transition matrix
#' Wcol <- gsemb_transition_matrix(adj, normalize = "col")
#' Matrix::colSums(Wcol)
#'
#' # Row-normalized transition matrix
#' Wrow <- gsemb_transition_matrix(adj, normalize = "row")
#' Matrix::rowSums(Wrow)
#' @export
gsemb_transition_matrix <- function(adj, normalize = c("col", "row"), eps = 1e-12) {
  normalize <- match.arg(normalize)
  if (!inherits(adj, "Matrix")) stop("adj must be a Matrix sparse/dense matrix")
  if (normalize == "col") {
    s <- Matrix::colSums(adj)
    s <- pmax(as.numeric(s), eps)
    Dinv <- Matrix::Diagonal(x = 1 / s)
    W <- adj %*% Dinv
  } else {
    s <- Matrix::rowSums(adj)
    s <- pmax(as.numeric(s), eps)
    Dinv <- Matrix::Diagonal(x = 1 / s)
    W <- Dinv %*% adj
  }
  Matrix::drop0(W)
}

#' Select landmark nodes for diffusion features
#'
#' Choose landmark nodes either by degree (highest first) or random sampling.
#'
#' @param adj Adjacency matrix with rownames as node IDs.
#' @param k Number of landmarks to select (capped at number of nodes).
#' @param method \code{"degree"} or \code{"random"}.
#' @param seed Random seed used when \code{method="random"}.
#'
#' @return A character vector of landmark node IDs.
#' @examples
#' # Build a small graph
#' edges <- data.frame(
#'   node1 = c("A", "B", "C", "D", "E"),
#'   node2 = c("B", "C", "D", "E", "A")
#' )
#' adj <- gsemb_build_graph(edges)
#'
#' # Select 2 landmarks by degree
#' lm_deg <- gsemb_select_landmarks(adj, k = 2, method = "degree")
#' lm_deg
#'
#' # Select 2 landmarks randomly
#' lm_rnd <- gsemb_select_landmarks(adj, k = 2, method = "random", seed = 42)
#' lm_rnd
#' @export
gsemb_select_landmarks <- function(adj, k = 128, method = c("degree", "random"), seed = 1) {
  method <- match.arg(method)
  nodes <- rownames(adj)
  if (is.null(nodes)) stop("adj must have rownames")
  n <- length(nodes)
  k <- min(k, n)
  if (method == "degree") {
    deg <- Matrix::rowSums(adj != 0)
    ord <- order(deg, decreasing = TRUE)
    nodes[ord][seq_len(k)]
  } else {
    set.seed(seed)
    sample(nodes, k)
  }
}
