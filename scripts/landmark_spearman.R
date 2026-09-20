# Source after full-node scores are built:
#   landmark_spearman(scores)
#
# scores: named list of named numeric vectors (one score per node).
# Methods without a global score can pass NULL.

landmark_spearman <- function(scores) {
  if (!is.list(scores) || !length(scores) || is.null(names(scores)) || any(names(scores) == "")) {
    stop("scores must be a non-empty named list", call. = FALSE)
  }

  named_scores <- Filter(Negate(is.null), scores)
  common_nodes <- if (length(named_scores)) {
    Reduce(intersect, lapply(named_scores, function(x) {
      if (is.null(names(x))) stop("each non-NULL score must be a named numeric vector", call. = FALSE)
      names(x)
    }))
  } else {
    character()
  }

  pairs <- utils::combn(names(scores), 2L, simplify = FALSE)
  spearman_tbl <- do.call(rbind, lapply(pairs, function(p) {
    x <- scores[[p[[1]]]]
    y <- scores[[p[[2]]]]
    reason <- NA_character_
    rho <- NA_real_
    n_nodes <- length(common_nodes)
    if (is.null(x) || is.null(y)) {
      reason <- "method_has_no_global_score"
    } else if (!length(common_nodes)) {
      reason <- "no_common_nodes"
    } else if (stats::sd(x[common_nodes]) == 0 || stats::sd(y[common_nodes]) == 0) {
      reason <- "constant_score"
    } else {
      rho <- suppressWarnings(stats::cor(x[common_nodes], y[common_nodes], method = "spearman"))
    }
    data.frame(
      method_a = p[[1]],
      method_b = p[[2]],
      n_nodes = n_nodes,
      spearman = rho,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }))
  rownames(spearman_tbl) <- NULL
  spearman_tbl
}
