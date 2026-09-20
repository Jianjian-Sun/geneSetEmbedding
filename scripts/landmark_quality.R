# Source after landmarks are selected:
#   landmark_quality(landmarks, adj_all, adj_high)
#
# landmarks: named list of character landmark IDs
# adj_all / adj_high: sparse adjacency matrices (e.g. graphs$all, graphs$high_conf)

landmark_quality <- function(landmarks, adj_all, adj_high) {
  if (!is.list(landmarks) || !length(landmarks) ||
      is.null(names(landmarks)) || any(names(landmarks) == "")) {
    stop("landmarks must be a non-empty named list", call. = FALSE)
  }

  nodes_all <- rownames(adj_all)
  if (is.null(nodes_all)) stop("adj_all must have rownames", call. = FALSE)

  degree_all <- setNames(as.numeric(Matrix::rowSums(adj_all != 0)), nodes_all)
  strength_all <- setNames(as.numeric(Matrix::rowSums(adj_all)), nodes_all)
  mean_incident_confidence <- ifelse(
    degree_all > 0,
    strength_all / degree_all,
    NA_real_
  )
  names(mean_incident_confidence) <- nodes_all

  degree_high <- setNames(numeric(length(nodes_all)), nodes_all)
  nodes_high <- rownames(adj_high)
  if (!is.null(nodes_high) && length(nodes_high)) {
    high_deg <- as.numeric(Matrix::rowSums(adj_high != 0))
    names(high_deg) <- nodes_high
    overlap <- intersect(nodes_all, nodes_high)
    degree_high[overlap] <- high_deg[overlap]
  }

  median_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else stats::median(x)
  }

  summary <- do.call(rbind, lapply(names(landmarks), function(method) {
    lm <- unique(as.character(landmarks[[method]]))
    missing <- setdiff(lm, nodes_all)
    if (length(missing)) {
      stop(method, " has landmarks not in adj_all: ",
           paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    }
    data.frame(
      method = method,
      n_landmarks = length(lm),
      median_total_degree = median_or_na(degree_all[lm]),
      median_high_degree = median_or_na(degree_high[lm]),
      median_incident_confidence = median_or_na(mean_incident_confidence[lm]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary
}
