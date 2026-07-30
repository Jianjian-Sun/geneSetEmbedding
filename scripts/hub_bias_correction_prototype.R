# =============================================================================
# Prototype retained for reference. Production API now lives in:
#   R/weighted_enrichment.R
#     gsemb_weighted_gsea(..., degree_beta=), gsemb_weighted_ora(...),
#     gsemb_sweep_hub_beta()
# Prefer library(geneSetEmbedding) / devtools::load_all(".") over sourcing this
# file. This script is kept for historical / offline experiments only.
# =============================================================================
#
# Hub gene bias correction prototype (standalone). Soft membership reweighting:
#   w' ∝ w / degree(g)^beta
# has been merged into R/weighted_enrichment.R.
#
# 申请书问题
# -----------
# Hub 基因在 RWR / soft membership 中可能获得偏高权重。
# 校正形式：对 soft membership 做
#   w'_gS ∝ w_gS / degree(g)^β
# 再按列重归一化；β = 0 表示不校正（与现逻辑一致）。
#
# 本脚本做什么
# -------------
# 1. source 同目录的 weighted_ora_gsea_prototype.R（复用 soft-W / 置换 / ORA·GSEA）
# 2. 在 W 算完之后、T_S 之前插入度数惩罚
# 3. 提供带 degree_beta 的 Weighted-ORA / Weighted-GSEA 入口
# 4. 提供 β ∈ {0, 0.25, 0.5, 1} 的扫参比较助手
#
# 不做
# ----
# 不修改 R/api.R、R/concise.R、R/graph.R；合回时对照下方地图。
#
# 用法（在包根目录）
# -----------------
#   devtools::load_all(".")
#   source("scripts/hub_bias_correction_prototype.R", encoding = "UTF-8")
#   # gsemb_weighted_gsea_hub(..., degree_beta = 0.5)
#   # gsemb_sweep_hub_beta(r, fit, betas = c(0, 0.25, 0.5, 1), ...)
#
# -----------------------------------------------------------------------------
# 合回原脚本修改地图
# -----------------------------------------------------------------------------
#
# | 本脚本符号                      | 合回位置                         | 应做的事 |
# |--------------------------------|----------------------------------|----------|
# | .proto_node_degree_from_adj    | 新建于 R/00-utils.R 或 api.R 内部 | 从 fit$adj 取命名度数向量 |
# | .proto_apply_degree_beta       | 紧接 soft membership W 之后      | W / deg^β + 列归一化 |
# | degree_beta 参数               | gsemb_embedding_enrichment 及    | 默认 0；写入 roxygen |
# |                                | 未来 gsemb_weighted_*            | |
# | gsemb_sweep_hub_beta           | 分析脚本或 vignette，非必须入包  | 扫参比较 |
# | 测试                           | tests/testthat/                  | β=0 与无校正一致；β>0 改变 ES |
#
# =============================================================================

# ---- 定位并加载 Weighted 原型 ------------------------------------------------

.hub_proto_find_weighted_script <- function() {
  candidates <- c(
    "scripts/weighted_ora_gsea_prototype.R",
    file.path(getwd(), "scripts", "weighted_ora_gsea_prototype.R"),
    file.path("e:/genesetembedding/geneSetEmbedding/scripts", "weighted_ora_gsea_prototype.R")
  )
  # 若本文件路径可知，优先同目录
  ofile <- NULL
  if (exists("sys.frames") && length(sys.frames()) >= 1L) {
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f) && nzchar(f)) {
        ofile <- f
        break
      }
    }
  }
  if (!is.null(ofile)) {
    candidates <- c(
      file.path(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)),
                "weighted_ora_gsea_prototype.R"),
      candidates
    )
  }
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  stop(
    "找不到 scripts/weighted_ora_gsea_prototype.R；",
    "请在包根目录 source 本脚本，或先手动 source Weighted 原型。"
  )
}

if (!exists(".proto_soft_membership_W", mode = "function")) {
  weighted_script <- .hub_proto_find_weighted_script()
  source(weighted_script, encoding = "UTF-8", local = FALSE)
}

if (!exists(".proto_soft_membership_W", mode = "function")) {
  stop("Weighted 原型未成功加载：缺少 .proto_soft_membership_W")
}

# ---- 度数与 hub 校正 ---------------------------------------------------------

#' 从邻接矩阵提取命名度数向量（与 genes 对齐）
#'
#' 无行和近似无向度数；缺失节点记为 0，随后用 pmax(deg, 1) 避免除零。
.proto_node_degree_from_adj <- function(adj, genes) {
  if (is.null(adj)) {
    stop("需要 fit$adj（gsemb_embedding 对象应含 adj）才能做 hub 校正")
  }
  if (is.null(rownames(adj))) {
    stop("adj 必须有 rownames（基因/节点名）")
  }
  # Matrix 或 dense
  if (inherits(adj, "Matrix")) {
    deg_all <- Matrix::rowSums(adj)
  } else {
    deg_all <- rowSums(as.matrix(adj))
  }
  names(deg_all) <- rownames(adj)
  deg <- as.numeric(deg_all[genes])
  names(deg) <- genes
  deg[is.na(deg)] <- 0
  deg
}

#' Soft membership hub 校正：W <- W / degree^β，再列归一化
#'
#' @param W n_genes × n_sets 权重矩阵（已 softmax）
#' @param degree 与 rownames(W) 同序的度数（可含 0）
#' @param beta 惩罚强度；0 = 不改动（数值上应与输入 W 一致，仍会走一遍归一化）
#' @param eps 数值稳定
.proto_apply_degree_beta <- function(W, degree, beta = 0, eps = 1e-8) {
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta) || beta < 0) {
    stop("degree_beta (beta) 必须是有限非负标量")
  }
  W <- as.matrix(W)
  genes <- rownames(W)
  if (is.null(genes)) stop("W 必须有 rownames")
  if (is.null(names(degree))) {
    if (length(degree) != nrow(W)) stop("degree 长度须等于 nrow(W)")
  } else {
    degree <- as.numeric(degree[genes])
    names(degree) <- genes
  }
  if (length(degree) != nrow(W)) stop("degree 与 W 行数不一致")

  if (isTRUE(all.equal(beta, 0))) {
    # β=0：不做惩罚；直接返回（避免无意义的重归一化浮点差）
    return(W)
  }

  deg_safe <- pmax(as.numeric(degree), 1)
  penalty <- deg_safe^beta
  W2 <- W / penalty
  W2[!is.finite(W2)] <- 0
  cs <- colSums(W2)
  W2 <- W2 / rep(pmax(cs, eps), each = nrow(W2))
  W2[, cs <= 0] <- 0
  W2[!is.finite(W2)] <- 0
  rownames(W2) <- genes
  colnames(W2) <- colnames(W)
  W2
}

#' Soft membership + 可选 hub 校正
.proto_soft_membership_W_hub <- function(gene_emb,
                                         mu,
                                         var,
                                         gene_sets = NULL,
                                         restrict_to_members = FALSE,
                                         score = c("loglik", "neg_mahalanobis"),
                                         temperature = 1.0,
                                         eps = 1e-8,
                                         degree = NULL,
                                         degree_beta = 0) {
  W <- .proto_soft_membership_W(
    gene_emb = gene_emb,
    mu = mu,
    var = var,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    eps = eps
  )
  if (is.null(degree)) {
    if (!isTRUE(all.equal(degree_beta, 0))) {
      stop("degree_beta > 0 时必须提供 degree 向量（或在驱动层传入 fit$adj）")
    }
    return(W)
  }
  .proto_apply_degree_beta(W, degree = degree, beta = degree_beta, eps = eps)
}

# ---- 带 hub 校正的单批检验 / 驱动 --------------------------------------------

.proto_weighted_test_sets_hub <- function(gene_stats,
                                          gene_emb,
                                          genes,
                                          mu,
                                          var,
                                          gene_sets,
                                          restrict_to_members,
                                          score,
                                          temperature,
                                          nperm,
                                          alternative,
                                          seed,
                                          eps,
                                          top_genes,
                                          perm_mat,
                                          degree,
                                          degree_beta) {
  W <- .proto_soft_membership_W_hub(
    gene_emb = gene_emb,
    mu = mu,
    var = var,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    eps = eps,
    degree = degree,
    degree_beta = degree_beta
  )

  stats_vec <- as.numeric(gene_stats)
  es <- as.numeric(crossprod(W, stats_vec))
  names(es) <- colnames(W)
  n_es <- length(es)

  z <- rep(NA_real_, n_es)
  pvals <- rep(NA_real_, n_es)

  if (nperm > 0L) {
    null_scores <- .proto_enrichment_null_scores(W, stats_vec, nperm, seed, perm_mat)
    null_mean <- colMeans(null_scores)
    null_sd <- pmax(apply(null_scores, 2, stats::sd), eps)
    z <- (es - null_mean) / null_sd

    es_mat <- matrix(es, nrow = nperm, ncol = n_es, byrow = TRUE)
    if (alternative == "greater") {
      pvals <- colMeans(null_scores >= es_mat)
    } else if (alternative == "less") {
      pvals <- colMeans(null_scores <= es_mat)
    } else {
      cen_mat <- matrix(null_mean, nrow = nperm, ncol = n_es, byrow = TRUE)
      shift_mat <- matrix(abs(es - null_mean), nrow = nperm, ncol = n_es, byrow = TRUE)
      pvals <- colMeans(abs(null_scores - cen_mat) >= shift_mat)
    }
    pvals <- pmax(pvals, 1 / nperm)
  }

  set_size <- rep(NA_integer_, n_es)
  if (!is.null(gene_sets)) {
    set_size <- vapply(names(es), function(sid) {
      if (!sid %in% names(gene_sets)) {
        return(NA_integer_)
      }
      length(intersect(gene_sets[[sid]], genes))
    }, integer(1))
  }

  core <- rep(NA_character_, length(es))
  if (top_genes > 0L) {
    for (j in seq_along(es)) {
      wj <- W[, j]
      ord <- order(wj, decreasing = TRUE)
      ord <- ord[seq_len(min(top_genes, length(ord)))]
      core[j] <- paste0(rownames(W)[ord], collapse = "/")
    }
  }

  data.frame(
    ID = names(es),
    ES = as.numeric(es),
    z = as.numeric(z),
    pvalue = as.numeric(pvals),
    p.adjust = NA_real_,
    setSize = as.integer(set_size),
    core_enrichment = core,
    degree_beta = as.numeric(degree_beta),
    stringsAsFactors = FALSE
  )
}

.proto_weighted_enrichment_driver_hub <- function(gene_stats,
                                                 x,
                                                 sets = NULL,
                                                 gene_sets = NULL,
                                                 restrict_to_members = FALSE,
                                                 score = c("loglik", "neg_mahalanobis"),
                                                 temperature = 1.0,
                                                 nperm = 1000,
                                                 alternative = c("two.sided", "greater", "less"),
                                                 seed = 1,
                                                 eps = 1e-8,
                                                 top_genes = 30,
                                                 degree_beta = 0,
                                                 adj = NULL) {
  score <- match.arg(score)
  alternative <- match.arg(alternative)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  if (!is.numeric(gene_stats) || is.null(names(gene_stats))) {
    stop("gene_stats must be a named numeric vector")
  }
  if (!is.numeric(degree_beta) || length(degree_beta) != 1L ||
        !is.finite(degree_beta) || degree_beta < 0) {
    stop("degree_beta must be a finite non-negative scalar")
  }
  if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 0) {
    stop("nperm must be a non-negative integer")
  }
  nperm <- as.integer(nperm)
  top_genes <- as.integer(top_genes)
  if (isTRUE(restrict_to_members) && is.null(gene_sets)) {
    stop("restrict_to_members=TRUE 时必须提供 gene_sets")
  }

  gene_emb <- x$gene_embedding
  if (is.null(rownames(gene_emb))) stop("gene embedding must have rownames")
  mu <- .proto_as_numeric_matrix(x$set_mu)
  var <- .proto_as_numeric_matrix(x$set_var)

  genes <- intersect(names(gene_stats), rownames(gene_emb))
  if (length(genes) < 2) {
    stop("not enough genes overlap between gene_stats and gene_embedding")
  }
  gene_emb <- gene_emb[genes, , drop = FALSE]
  gene_stats <- gene_stats[rownames(gene_emb)]
  storage.mode(gene_stats) <- "double"
  genes <- rownames(gene_emb)

  # 度数：优先显式 adj，否则 fit$adj
  if (is.null(adj)) {
    adj <- x$adj
  }
  degree <- NULL
  if (!isTRUE(all.equal(degree_beta, 0))) {
    degree <- .proto_node_degree_from_adj(adj, genes)
  } else if (!is.null(adj)) {
    # β=0 也可预计算，便于与扫参共用；非必须
    degree <- .proto_node_degree_from_adj(adj, genes)
  }

  set_names <- rownames(mu)
  if (is.null(sets)) {
    set_idx <- seq_len(nrow(mu))
  } else {
    set_idx <- match(as.character(sets), set_names, nomatch = 0L)
    if (any(set_idx == 0L)) {
      bad <- sets[set_idx == 0L]
      stop(
        "some sets are not in set_mu rownames; first missing: ",
        paste(utils::head(bad, 3L), collapse = ", ")
      )
    }
  }
  if (length(set_idx) == 0L) stop("no sets found in embedding")

  mu <- mu[set_idx, , drop = FALSE]
  var <- var[set_idx, , drop = FALSE]
  n_local <- nrow(mu)
  local_idx <- seq_len(n_local)

  if (!is.null(gene_sets)) {
    gene_sets <- .proto_validate_gene_sets(gene_sets)
  }

  perm_mat <- NULL
  if (nperm > 0L) {
    stats_vec <- as.numeric(gene_stats)
    set.seed(seed)
    perm_mat <- vapply(
      seq_len(nperm),
      function(b) sample(stats_vec, length(stats_vec), replace = FALSE),
      FUN.VALUE = numeric(length(stats_vec))
    )
  }

  n_genes <- nrow(gene_emb)
  sets_per_batch <- if (as.numeric(n_genes) * as.numeric(n_local) > 2e7) {
    max(50L, as.integer(floor(2e7 / max(n_genes, 1L))))
  } else {
    n_local
  }
  if (!is.finite(sets_per_batch) || sets_per_batch < 1L) {
    sets_per_batch <- n_local
  }

  idx_batches <- split(local_idx, ceiling(seq_along(local_idx) / sets_per_batch))
  parts <- lapply(idx_batches, function(batch_idx) {
    .proto_weighted_test_sets_hub(
      gene_stats = gene_stats,
      gene_emb = gene_emb,
      genes = genes,
      mu = mu[batch_idx, , drop = FALSE],
      var = var[batch_idx, , drop = FALSE],
      gene_sets = gene_sets,
      restrict_to_members = restrict_to_members,
      score = score,
      temperature = temperature,
      nperm = nperm,
      alternative = alternative,
      seed = seed,
      eps = eps,
      top_genes = top_genes,
      perm_mat = perm_mat,
      degree = degree,
      degree_beta = degree_beta
    )
  })

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out$p.adjust <- stats::p.adjust(out$pvalue, method = "BH")
  out
}

# ---- 对外 API ----------------------------------------------------------------

#' Weighted-ORA + hub 校正（原型）
#'
#' @param degree_beta 度数惩罚 β；0 = 不校正。常用 0, 0.25, 0.5, 1。
#' @param adj 可选邻接矩阵；默认用 x$adj。
gsemb_weighted_ora_hub <- function(y,
                                   x,
                                   sets = NULL,
                                   gene_sets = NULL,
                                   restrict_to_members = FALSE,
                                   degree_beta = 0,
                                   adj = NULL,
                                   score = c("loglik", "neg_mahalanobis"),
                                   temperature = 1.0,
                                   nperm = 1000,
                                   alternative = c("two.sided", "greater", "less"),
                                   seed = 1,
                                   eps = 1e-8,
                                   top_genes = 30) {
  if (!is.numeric(y) || is.null(names(y))) {
    stop("y must be a named numeric vector (0/1 labels)")
  }
  y_names <- names(y)
  y_num <- as.numeric(y)
  names(y_num) <- y_names
  uniq <- unique(y_num[is.finite(y_num)])
  if (length(uniq) == 0L || !all(uniq %in% c(0, 1))) {
    stop("Weighted-ORA 要求 y 仅含 0 和 1")
  }
  storage.mode(y_num) <- "double"

  .proto_weighted_enrichment_driver_hub(
    gene_stats = y_num,
    x = x,
    sets = sets,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    nperm = nperm,
    alternative = alternative,
    seed = seed,
    eps = eps,
    top_genes = top_genes,
    degree_beta = degree_beta,
    adj = adj
  )
}

#' Weighted-GSEA + hub 校正（原型）
gsemb_weighted_gsea_hub <- function(r,
                                    x,
                                    sets = NULL,
                                    gene_sets = NULL,
                                    restrict_to_members = FALSE,
                                    degree_beta = 0,
                                    adj = NULL,
                                    score = c("loglik", "neg_mahalanobis"),
                                    temperature = 1.0,
                                    nperm = 1000,
                                    alternative = c("two.sided", "greater", "less"),
                                    seed = 1,
                                    eps = 1e-8,
                                    top_genes = 30) {
  if (!is.numeric(r) || is.null(names(r))) {
    stop("r must be a named numeric vector")
  }

  .proto_weighted_enrichment_driver_hub(
    gene_stats = r,
    x = x,
    sets = sets,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    nperm = nperm,
    alternative = alternative,
    seed = seed,
    eps = eps,
    top_genes = top_genes,
    degree_beta = degree_beta,
    adj = adj
  )
}

#' 扫多个 β，比较校正前后富集（便于申请书实验）
#'
#' @param gene_stats 命名向量；mode="ora" 时须为 0/1，mode="gsea" 时为连续。
#' @param betas 数值向量，默认 c(0, 0.25, 0.5, 1)。
#' @param mode "gsea" 或 "ora"。
#' @return list(results = 按 β 的 data.frame 列表, combined = 带 degree_beta 列的长表,
#'   compare = 各 set 在不同 β 下的 ES 宽表)。
gsemb_sweep_hub_beta <- function(gene_stats,
                                 x,
                                 betas = c(0, 0.25, 0.5, 1),
                                 mode = c("gsea", "ora"),
                                 sets = NULL,
                                 gene_sets = NULL,
                                 restrict_to_members = FALSE,
                                 adj = NULL,
                                 nperm = 200,
                                 seed = 1,
                                 alternative = "greater",
                                 ...) {
  mode <- match.arg(mode)
  betas <- as.numeric(betas)
  if (any(!is.finite(betas) | betas < 0)) stop("betas 必须全为有限非负数")

  fun <- if (mode == "ora") gsemb_weighted_ora_hub else gsemb_weighted_gsea_hub
  results <- lapply(betas, function(b) {
    fun(
      gene_stats, x,
      sets = sets,
      gene_sets = gene_sets,
      restrict_to_members = restrict_to_members,
      degree_beta = b,
      adj = adj,
      nperm = nperm,
      seed = seed,
      alternative = alternative,
      ...
    )
  })
  names(results) <- paste0("beta_", betas)

  combined <- do.call(rbind, results)
  rownames(combined) <- NULL

  # ES 宽表：行 = set，列 = beta_*
  es_list <- lapply(seq_along(betas), function(i) {
    d <- results[[i]]
    v <- d$ES
    names(v) <- d$ID
    v
  })
  all_ids <- unique(unlist(lapply(es_list, names), use.names = FALSE))
  compare <- data.frame(ID = all_ids, stringsAsFactors = FALSE)
  for (i in seq_along(betas)) {
    compare[[paste0("ES_beta_", betas[i])]] <- as.numeric(es_list[[i]][all_ids])
  }
  if ("ES_beta_0" %in% names(compare)) {
    for (b in betas[betas != 0]) {
      col <- paste0("ES_beta_", b)
      compare[[paste0("deltaES_vs0_", b)]] <- compare[[col]] - compare$ES_beta_0
    }
  }

  list(results = results, combined = combined, compare = compare)
}

# ---- 玩具演示（默认不执行）--------------------------------------------------

if (FALSE) {
  # devtools::load_all(".")
  # source("scripts/hub_bias_correction_prototype.R", encoding = "UTF-8")

  set.seed(1)
  edges <- data.frame(
    node1 = c("A", "A", "A", "B", "C", "D"),
    node2 = c("B", "C", "D", "C", "D", "E"),
    weight = c(1, 1, 1, 1, 1, 1)
  )
  # A 是 hub（度数更高）
  gene_sets <- list(S1 = c("A", "B", "C"), S2 = c("C", "D", "E"))
  fit <- gsemb_fit(
    edges, gene_sets,
    weight = "weight", method = "svd",
    dim = 2, k = 2, max_iter = 30, epochs = 3
  )

  r <- c(A = 2, B = 1, C = 0.5, D = -1, E = -2)
  sweep <- gsemb_sweep_hub_beta(
    r, fit,
    betas = c(0, 0.25, 0.5, 1),
    mode = "gsea",
    gene_sets = gene_sets,
    nperm = 50, seed = 1,
    alternative = "greater"
  )
  print(sweep$compare)

  # β=0 应与无 hub 的 Weighted-GSEA 原型一致
  g0 <- gsemb_weighted_gsea_hub(r, fit, gene_sets = gene_sets, degree_beta = 0, nperm = 50, seed = 1, alternative = "greater")
  g_ref <- gsemb_weighted_gsea(r, fit, gene_sets = gene_sets, nperm = 50, seed = 1, alternative = "greater")
  cat("beta0 vs weighted prototype ES diff:", max(abs(g0$ES - g_ref$ES)), "\n")
}

message(
  "已加载 Hub bias 原型：gsemb_weighted_ora_hub() / gsemb_weighted_gsea_hub() / ",
  "gsemb_sweep_hub_beta()；原 R/ 未改动。"
)
