# Source after baseline landmarks are selected:
#   landmark_dropout(landmarks, edges, k = 128L)
#
# Randomly drops dropout_fraction of undirected edges, rebuilds graphs,
# reselects landmarks with the same four method rules, and reports Jaccard
# vs the baseline sets.

landmark_dropout <- function(landmarks,
                             edges,
                             k,
                             node1 = "protein1",
                             node2 = "protein2",
                             weight = "combined_score",
                             high_cut = 700,
                             dropout_fraction = 0.05,
                             repeats = 100L,
                             seed = 20260917L,
                             selectors = NULL) {
  if (!is.list(landmarks) || !length(landmarks) ||
      is.null(names(landmarks)) || any(names(landmarks) == "")) {
    stop("landmarks must be a non-empty named list", call. = FALSE)
  }
  if (is.null(selectors)) {
    selectors <- list(
      degree_all = list(graph = "all", method = "degree"),
      weighted_degree_all = list(graph = "all", method = "weighted_degree"),
      degree_high = list(graph = "high", method = "degree"),
      weighted_degree_high = list(graph = "high", method = "weighted_degree")
    )
  }
  selectors <- selectors[intersect(names(selectors), names(landmarks))]
  if (!length(selectors)) stop("no selectors match landmarks names", call. = FALSE)

  jaccard <- function(a, b) {
    a <- unique(a)
    b <- unique(b)
    u <- union(a, b)
    if (!length(u)) return(NA_real_)
    length(intersect(a, b)) / length(u)
  }

  a <- trimws(as.character(edges[[node1]]))
  b <- trimws(as.character(edges[[node2]]))
  w <- as.numeric(edges[[weight]])
  ok <- !is.na(a) & !is.na(b) & a != "" & b != "" & a != b & is.finite(w) & w >= 0
  left <- pmin(a[ok], b[ok])
  right <- pmax(a[ok], b[ok])
  w <- w[ok]
  pair <- paste(left, right, sep = "\r")
  w <- ave(w, pair, FUN = max)
  keep <- !duplicated(pair)
  clean <- data.frame(
    node1 = left[keep], node2 = right[keep], weight = w[keep],
    stringsAsFactors = FALSE
  )
  n_edges <- nrow(clean)
  n_drop <- if (dropout_fraction <= 0 || n_edges == 0) {
    0L
  } else {
    max(1L, floor(n_edges * dropout_fraction))
  }

  select_on_graphs <- function(adj_all, adj_high) {
    lapply(selectors, function(sel) {
      adj <- if (identical(sel$graph, "high")) adj_high else adj_all
      gsemb_select_landmarks(adj, k = k, method = sel$method)
    })
  }

  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, repeats)
  rows <- vector("list", repeats * length(selectors))
  at <- 1L
  for (replicate_id in seq_len(repeats)) {
    set.seed(seeds[[replicate_id]])
    dropped <- if (n_drop) sort(sample.int(n_edges, n_drop)) else integer()
    kept <- if (length(dropped)) clean[-dropped, , drop = FALSE] else clean
    high <- kept[kept$weight > high_cut, , drop = FALSE]
    adj_all <- gsemb_build_graph(
      kept, node1 = "node1", node2 = "node2", weight = "weight", directed = FALSE
    )
    adj_high <- gsemb_build_graph(
      high, node1 = "node1", node2 = "node2", weight = "weight", directed = FALSE
    )
    selected <- select_on_graphs(adj_all, adj_high)
    for (method in names(selectors)) {
      rows[[at]] <- data.frame(
        replicate = replicate_id,
        seed = seeds[[replicate_id]],
        method = method,
        n_edges_removed = length(dropped),
        overlap_n = length(intersect(landmarks[[method]], selected[[method]])),
        jaccard = jaccard(landmarks[[method]], selected[[method]]),
        stringsAsFactors = FALSE
      )
      at <- at + 1L
    }
  }
  replicates <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(replicates$jaccard, replicates$method), function(x) {
    data.frame(
      mean = mean(x),
      sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
      median = stats::median(x),
      q025 = unname(stats::quantile(x, 0.025)),
      q975 = unname(stats::quantile(x, 0.975)),
      min = min(x),
      max = max(x)
    )
  }))
  summary$method <- rownames(summary)
  summary <- summary[, c("method", setdiff(names(summary), "method")), drop = FALSE]
  rownames(summary) <- NULL
  list(replicates = replicates, summary = summary)
}
