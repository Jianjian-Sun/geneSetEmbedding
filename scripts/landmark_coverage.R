# Source after landmarks are selected:
#   landmark_coverage(landmarks, adj_eval, radius = 2L)
#
# landmarks: named list of character landmark IDs
# adj_eval: adjacency used for hop distances (usually graphs$high_conf)

landmark_coverage <- function(landmarks, adj_eval, radius = 2L) {
  if (!is.list(landmarks) || !length(landmarks) ||
      is.null(names(landmarks)) || any(names(landmarks) == "")) {
    stop("landmarks must be a non-empty named list", call. = FALSE)
  }
  nodes <- rownames(adj_eval)
  if (is.null(nodes)) stop("adj_eval must have rownames", call. = FALSE)
  radius <- as.integer(radius)

  sm <- Matrix::summary(adj_eval)
  neighbors <- setNames(vector("list", length(nodes)), nodes)
  if (nrow(sm)) {
    src <- nodes[sm$i]
    dst <- nodes[sm$j]
    split_n <- split(dst, src)
    neighbors[names(split_n)] <- lapply(split_n, unique)
  }

  nearest_landmark_distance <- function(landmarks) {
    distance <- setNames(rep(Inf, length(nodes)), nodes)
    seeds <- intersect(unique(as.character(landmarks)), nodes)
    if (!length(seeds)) return(distance)
    distance[seeds] <- 0
    queue <- seeds
    head <- 1L
    while (head <= length(queue)) {
      current <- queue[[head]]
      head <- head + 1L
      nbrs <- neighbors[[current]]
      unseen <- nbrs[is.infinite(distance[nbrs])]
      if (length(unseen)) {
        distance[unseen] <- distance[[current]] + 1
        queue <- c(queue, unseen)
      }
    }
    distance
  }

  summary <- do.call(rbind, lapply(names(landmarks), function(method) {
    d <- nearest_landmark_distance(landmarks[[method]])
    reachable <- d[is.finite(d)]
    data.frame(
      method = method,
      radius = radius,
      coverage_at_2 = mean(d <= 2L),
      coverage_at_radius = mean(d <= radius),
      mean_nearest_distance = if (length(reachable)) mean(reachable) else NA_real_,
      unreachable_fraction = mean(!is.finite(d)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary
}
