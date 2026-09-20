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
#' Choose landmark nodes by degree, weighted degree (node strength), random
#' sampling, betweenness centrality (unweighted or weight-aware), or K-Medoids
#' (hop distance or \code{1/weight} path costs; farthest-first then restricted PAM).
#'
#' @param adj Adjacency matrix with rownames as node IDs.
#' @param k Number of landmarks to select (capped at number of nodes).
#' @param method Landmark selection method: \code{"degree"} counts incident
#'   nonzero edges, \code{"weighted_degree"} sums incident edge weights,
#'   \code{"random"} samples nodes uniformly, \code{"betweenness"} ranks by
#'   unweighted shortest-path betweenness on \code{adj != 0},
#'   \code{"weighted_betweenness"} ranks by betweenness where edge lengths are
#'   \code{1 / weight} (high confidence / weight preferred as shorter),
#'   \code{"kmedoids"} selects hop-distance medoids on the largest connected
#'   component of \code{adj != 0}, and \code{"weighted_kmedoids"} uses the same
#'   procedure with edge costs \code{1 / weight}.
#' @param seed Random seed used when \code{method="random"} or either
#'   K-Medoids method (farthest-first start).
#' @param betweenness_cutoff Path-length cutoff for
#'   \code{igraph::betweenness} when \code{method} is \code{"betweenness"} or
#'   \code{"weighted_betweenness"}. Default \code{4}. For unweighted graphs this
#'   is hop length; for weighted it is on the sum of edge costs
#'   (\code{1/weight}). Use \code{-1} for exact (untruncated) betweenness.
#'   Ignored for other methods.
#' @param kmedoids_m Per-cluster candidate cap (degree Top-M for
#'   \code{"kmedoids"}, strength Top-M for \code{"weighted_kmedoids"}) for
#'   restricted PAM. Ignored for other methods.
#' @param kmedoids_max_iter Maximum PAM passes for either K-Medoids method.
#'   Ignored for other methods.
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
#' # Select 2 landmarks by weighted degree (node strength)
#' lm_weighted <- gsemb_select_landmarks(adj, k = 2, method = "weighted_degree")
#' lm_weighted
#'
#' # Select 2 landmarks randomly
#' lm_rnd <- gsemb_select_landmarks(adj, k = 2, method = "random", seed = 42)
#' lm_rnd
#'
#' # Select by truncated betweenness (default cutoff = 4)
#' lm_bc <- gsemb_select_landmarks(adj, k = 2, method = "betweenness")
#' lm_bc
#'
#' # Weight-aware betweenness (edge cost = 1/weight)
#' lm_wbc <- gsemb_select_landmarks(adj, k = 2, method = "weighted_betweenness")
#' lm_wbc
#'
#' # Hop-distance K-Medoids on a path (k = 1 converges to the center)
#' path <- data.frame(
#'   node1 = c("A", "B", "C", "D"),
#'   node2 = c("B", "C", "D", "E")
#' )
#' gsemb_select_landmarks(gsemb_build_graph(path), k = 1, method = "kmedoids")
#'
#' # Weight-aware K-Medoids (edge cost = 1/weight)
#' wpath <- data.frame(
#'   node1 = c("A", "B", "A", "D"),
#'   node2 = c("B", "C", "D", "C"),
#'   weight = c(100, 100, 1, 1)
#' )
#' gsemb_select_landmarks(
#'   gsemb_build_graph(wpath, weight = "weight"),
#'   k = 1,
#'   method = "weighted_kmedoids"
#' )
#' @export
gsemb_select_landmarks <- function(adj,
                                   k = 128,
                                   method = c(
                                     "degree",
                                     "weighted_degree",
                                     "random",
                                     "betweenness",
                                     "weighted_betweenness",
                                     "kmedoids",
                                     "weighted_kmedoids"
                                   ),
                                   seed = 1,
                                   betweenness_cutoff = 4,
                                   kmedoids_m = 50,
                                   kmedoids_max_iter = 10) {
  method <- match.arg(method)
  nodes <- rownames(adj)
  if (is.null(nodes)) stop("adj must have rownames")
  n <- length(nodes)
  k <- min(k, n)
  if (method == "degree") {
    deg <- Matrix::rowSums(adj != 0)
    ord <- order(deg, decreasing = TRUE)
    nodes[ord][seq_len(k)]
  } else if (method == "weighted_degree") {
    strength <- Matrix::rowSums(adj)
    ord <- order(strength, decreasing = TRUE)
    nodes[ord][seq_len(k)]
  } else if (method == "random") {
    set.seed(seed)
    sample(nodes, k)
  } else if (method == "betweenness") {
    g <- igraph::graph_from_adjacency_matrix(
      adj != 0,
      mode = "undirected",
      diag = FALSE
    )
    bc <- igraph::betweenness(g, directed = FALSE, cutoff = betweenness_cutoff)
    ord <- order(bc, decreasing = TRUE)
    nodes[ord][seq_len(k)]
  } else if (method == "kmedoids" || method == "weighted_kmedoids") {
    .select_landmarks_kmedoids(
      adj,
      k = k,
      seed = seed,
      m = kmedoids_m,
      max_iter = kmedoids_max_iter,
      weighted = identical(method, "weighted_kmedoids")
    )
  } else {
    # ponytail: igraph treats weights as path costs; invert so high score = short
    g <- igraph::graph_from_adjacency_matrix(
      adj,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    costs <- 1 / pmax(igraph::E(g)$weight, .Machine$double.eps)
    bc <- igraph::betweenness(
      g,
      directed = FALSE,
      weights = costs,
      cutoff = betweenness_cutoff
    )
    ord <- order(bc, decreasing = TRUE)
    nodes[ord][seq_len(k)]
  }
}

# ponytail: largest CC only; split K across components if small CCs must be covered.
.select_landmarks_kmedoids <- function(adj, k, seed, m, max_iter, weighted = FALSE) {
  g <- igraph::graph_from_adjacency_matrix(
    adj != 0,
    mode = "undirected",
    diag = FALSE
  )
  comp <- igraph::components(g)
  keep <- which(comp$membership == which.max(comp$csize))
  nodes <- igraph::V(g)$name[keep]
  n <- length(nodes)
  k <- min(as.integer(k), n)
  if (k >= n) {
    return(nodes)
  }

  # CSC of the induced subgraph for C++ SSSP / restricted PAM.
  sub <- Matrix::drop0(adj[nodes, nodes, drop = FALSE])
  idx0 <- kmedoids_hop(
    as.integer(sub@i),
    as.integer(sub@p),
    as.numeric(sub@x),
    as.integer(k),
    as.integer(seed),
    as.integer(m),
    as.integer(max_iter),
    isTRUE(weighted)
  )
  nodes[as.integer(idx0) + 1L]
}
