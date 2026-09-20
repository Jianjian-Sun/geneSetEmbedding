# Source after landmarks are selected:
#   landmark_jaccard(landmarks)

landmark_jaccard <- function(landmarks) {
  jaccard <- function(a, b) {
    a <- unique(a)
    b <- unique(b)
    u <- union(a, b)
    if (!length(u)) return(NA_real_)
    length(intersect(a, b)) / length(u)
  }

  pairs <- utils::combn(names(landmarks), 2L, simplify = FALSE)
  jaccard_tbl <- do.call(rbind, lapply(pairs, function(p) {
    a <- landmarks[[p[[1]]]]
    b <- landmarks[[p[[2]]]]
    data.frame(
      method_a = p[[1]],
      method_b = p[[2]],
      intersection_n = length(intersect(a, b)),
      union_n = length(union(a, b)),
      jaccard = jaccard(a, b)
    )
  }))
  rownames(jaccard_tbl) <- NULL
  jaccard_tbl
}
