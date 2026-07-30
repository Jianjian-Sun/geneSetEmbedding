# =============================================================================
# Prototype retained for reference. Production API now lives in:
#   R/weighted_enrichment.R
#     gsemb_embedding_enrichment(), gsemb_weighted_gsea(), gsemb_weighted_ora()
# Prefer library(geneSetEmbedding) / devtools::load_all(".") over sourcing this
# file. This script is kept for historical / offline experiments only.
# =============================================================================
#
# Weighted-ORA / Weighted-GSEA (standalone prototype; logic merged into the package)
#
# 目的
# -----
# 申请书要求 soft membership 加权检验，并拆成：
#   - Weighted-ORA ：输入二值标签 y_g，T_S = sum_g w_gS * y_g
#   - Weighted-GSEA：输入连续统计量 r_g（t / logFC 等），T_S = sum_g w_gS * r_g
# 现有包函数 gsemb_embedding_enrichment()（现位于 R/weighted_enrichment.R）已具备
# 「连续统计量 + soft weight + 置换 + BH」；正式拆分与 hub 校正已合入包。
#
# 本文件是早期「可运行副本」；新工作请直接用包导出函数。
#
# 依赖
# -----
# 复用已安装/已 load_all 的 geneSetEmbedding：
#   - gsemb_gene_to_set_score（导出）
#   - as_numeric_matrix / validate_gene_sets（内部，通过 ::: 调用）
# 不在此重复实现打分，避免与 R/concise.R 双份算法。
#
# 用法
# -----
#   # 在包根目录：
#   devtools::load_all(".")
#   source("scripts/weighted_ora_gsea_prototype.R", encoding = "UTF-8")
#   # 然后调用 gsemb_weighted_ora() / gsemb_weighted_gsea()
#   # 文末 if (FALSE) { ... } 为玩具演示，默认不执行。
#
# -----------------------------------------------------------------------------
# 合回原脚本修改地图（确认后合入时对照此表；本次不改这些文件）
# -----------------------------------------------------------------------------
#
# | 副本中的符号                    | 对应原文件位置                                      | 合回时应做的事 |
# |--------------------------------|-----------------------------------------------------|----------------|
# | .proto_soft_membership_W       | R/api.R  .gsemb_embedding_enrichment_sets           | 抽成共享内部函数；约 521–535 行打分→softmax→W；ORA/GSEA 共用 |
# |                                | （打分→softmax→W）                                  | 并在此接入 restrict_to_members |
# | .proto_enrichment_null_scores  | R/api.R  .gsemb_enrichment_null_scores 约 471–502 行 | 可原样复用（不必改名）；副本仅避免与包符号冲突 |
# | .proto_weighted_test_sets      | R/api.R  .gsemb_embedding_enrichment_sets 约 507–598 | 改为调用 soft-W + 置换；gene_sets 注释改为「也可用于 restrict」 |
# | gsemb_weighted_gsea            | R/api.R  gsemb_embedding_enrichment 约 649–762 行   | 重命名或薄封装；gene_stats 连续；新增 restrict_to_members |
# | gsemb_weighted_ora             | （原包无对应 API）                                   | 新增；校验 y∈{0,1}；内部调同一内核 |
# | restrict_to_members 分支       | 现逻辑无                                            | 建 W 时限制 softmax 支撑集到成员并归一化；更新 roxygen |
# | 导出                           | NAMESPACE                                           | export(gsemb_weighted_ora); export(gsemb_weighted_gsea)； |
# |                                |                                                     | 可选保留 gsemb_embedding_enrichment 为 GSEA 别名 |
# | 测试                           | tests/testthat/test-basic.R                         | 增加二值 ORA、连续 GSEA、restrict TRUE/FALSE 对照 |
#
# 兼容约定（自检）
# ---------------
# gsemb_weighted_gsea(..., restrict_to_members = FALSE) 在相同 seed/nperm/
# score/temperature/alternative 下，ES / pvalue 应与现
# gsemb_embedding_enrichment(...) 数值一致（允许浮点误差）。
#
# restrict_to_members
# -------------------
# - FALSE（默认）：与当前 gsemb_embedding_enrichment 一致。
#   支撑集 G_S = 全部「有统计量且有嵌入」的基因；非成员也参与 soft weight 与 T_S。
# - TRUE：仅 gene_sets[[S]] 与嵌入/统计量的交集进入该列 softmax 归一化；
#   非成员权重为 0。此时必须提供 gene_sets。
#
# 统计定义
# --------
#   w_gS = exp(s_gS / T) / sum_{g' in G_S} exp(s_g'S / T)
#   T_S  = sum_g w_gS * stats_g
#   置换 stats → 经验 p → 全集合并后 BH-FDR。
# =============================================================================

# ---- 依赖检查 ----------------------------------------------------------------

if (!requireNamespace("geneSetEmbedding", quietly = TRUE) &&
      !exists("gsemb_gene_to_set_score", mode = "function")) {
  stop(
    "请先加载 geneSetEmbedding：devtools::load_all('.') 或 library(geneSetEmbedding)"
  )
}

# 内部工具：优先用命名空间，兼容 load_all 后已挂到搜索路径的情况
.proto_as_numeric_matrix <- function(x, id_col = NULL) {
  if (exists("as_numeric_matrix", mode = "function", inherits = TRUE)) {
    return(as_numeric_matrix(x, id_col = id_col))
  }
  geneSetEmbedding:::as_numeric_matrix(x, id_col = id_col)
}

.proto_validate_gene_sets <- function(gene_sets) {
  if (exists("validate_gene_sets", mode = "function", inherits = TRUE)) {
    return(validate_gene_sets(gene_sets))
  }
  geneSetEmbedding:::validate_gene_sets(gene_sets)
}

.proto_gene_to_set_score <- function(...) {
  if (exists("gsemb_gene_to_set_score", mode = "function", inherits = TRUE)) {
    return(gsemb_gene_to_set_score(...))
  }
  geneSetEmbedding::gsemb_gene_to_set_score(...)
}

# ---- 副本：置换零分布（对应 api.R .gsemb_enrichment_null_scores）--------------

# 【副本说明】逻辑对齐原函数：优先用外层预生成的 perm_mat + crossprod；
# 否则尝试包内 Rcpp；再否则纯 R sample。不修改原 api.R。
.proto_enrichment_null_scores <- function(W,
                                          stats_vec,
                                          nperm,
                                          seed,
                                          perm_mat = NULL,
                                          n_workers = 1L) {
  stats_vec <- as.numeric(stats_vec)

  if (!is.null(perm_mat)) {
    return(t(crossprod(W, perm_mat)))
  }

  # 尝试包内已注册的 Rcpp 例程（若包已加载）
  if (isNamespaceLoaded("geneSetEmbedding")) {
    ns <- asNamespace("geneSetEmbedding")
    if (exists(".native_routine_available", envir = ns, inherits = FALSE) &&
          exists("rcpp_enrichment_permutations", envir = ns, inherits = FALSE)) {
      avail <- get(".native_routine_available", envir = ns)
      if (isTRUE(avail("_geneSetEmbedding_rcpp_enrichment_permutations"))) {
        return(get("rcpp_enrichment_permutations", envir = ns)(W, stats_vec, nperm, seed))
      }
    }
  }

  set.seed(seed)
  perm_inner <- vapply(
    seq_len(nperm),
    function(b) sample(stats_vec, length(stats_vec), replace = FALSE),
    FUN.VALUE = numeric(length(stats_vec))
  )

  if (n_workers <= 1L || !requireNamespace("future.apply", quietly = TRUE)) {
    return(t(crossprod(W, perm_inner)))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = n_workers)
  null_list <- future.apply::future_lapply(
    seq_len(nperm),
    function(b) as.numeric(crossprod(W, perm_inner[, b])),
    future.seed = TRUE
  )
  do.call(rbind, null_list)
}

# ---- 副本：soft membership 权重矩阵 W ----------------------------------------

# 对应合回：从 .gsemb_embedding_enrichment_sets 抽出的共享内核。
# restrict_to_members = TRUE 时，非成员分数置为 -Inf，再 softmax，
# 等价于仅在成员集合 G_S 上归一化（与申请书公式一致）。
.proto_soft_membership_W <- function(gene_emb,
                                     mu,
                                     var,
                                     gene_sets = NULL,
                                     restrict_to_members = FALSE,
                                     score = c("loglik", "neg_mahalanobis"),
                                     temperature = 1.0,
                                     eps = 1e-8) {
  score <- match.arg(score)
  if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) {
    stop("temperature must be a positive scalar")
  }
  if (isTRUE(restrict_to_members)) {
    if (is.null(gene_sets)) {
      stop("restrict_to_members=TRUE 时必须提供 gene_sets")
    }
    gene_sets <- .proto_validate_gene_sets(gene_sets)
  }

  S <- .proto_gene_to_set_score(
    gene_embedding = gene_emb,
    set_mu = mu,
    set_var = var,
    score = score,
    eps = eps
  )
  S <- as.matrix(S)
  genes <- rownames(gene_emb)
  set_ids <- rownames(mu)

  if (isTRUE(restrict_to_members)) {
    # 非成员不参与该 set 的 soft membership：softmax 前遮掉
    for (j in seq_along(set_ids)) {
      sid <- set_ids[j]
      if (!sid %in% names(gene_sets)) {
        S[, j] <- -Inf
        next
      }
      members <- intersect(gene_sets[[sid]], genes)
      if (length(members) == 0L) {
        S[, j] <- -Inf
        next
      }
      keep <- genes %in% members
      S[!keep, j] <- -Inf
    }
  }

  # 数值稳定的列 softmax
  col_max <- apply(S, 2, function(col) {
    finite <- col[is.finite(col)]
    if (length(finite) == 0L) return(0)
    max(finite)
  })
  S <- S - rep(col_max, each = nrow(S))
  exp_S <- exp(S / temperature)
  exp_S[!is.finite(exp_S)] <- 0
  col_sum <- colSums(exp_S)
  # 全被遮掉的列：权重全 0（后续 ES=0，p 需谨慎解读）
  W <- exp_S / rep(pmax(col_sum, eps), each = nrow(exp_S))
  W[, col_sum <= 0] <- 0
  W[!is.finite(W)] <- 0
  rownames(W) <- genes
  colnames(W) <- set_ids
  W
}

# ---- 副本：单批 sets 上的加权检验（ES + 置换 p）------------------------------

.proto_weighted_test_sets <- function(gene_stats,
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
                                      perm_mat) {
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
    stringsAsFactors = FALSE
  )
}

# ---- 共享驱动：对齐基因、选 sets、分批、BH -----------------------------------

.proto_weighted_enrichment_driver <- function(gene_stats,
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
                                              top_genes = 30) {
  score <- match.arg(score)
  alternative <- match.arg(alternative)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  if (!is.numeric(gene_stats) || is.null(names(gene_stats))) {
    stop("gene_stats must be a named numeric vector")
  }
  if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) {
    stop("temperature must be a positive scalar")
  }
  if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 0) {
    stop("nperm must be a non-negative integer")
  }
  nperm <- as.integer(nperm)
  if (!is.numeric(top_genes) || length(top_genes) != 1 || top_genes < 0) {
    stop("top_genes must be a non-negative integer")
  }
  top_genes <- as.integer(top_genes)
  if (isTRUE(restrict_to_members) && is.null(gene_sets)) {
    stop("restrict_to_members=TRUE 时必须提供 gene_sets")
  }

  gene_emb <- x$gene_embedding
  if (is.null(rownames(gene_emb))) stop("gene embedding must have rownames")
  mu <- .proto_as_numeric_matrix(x$set_mu)
  var <- .proto_as_numeric_matrix(x$set_var)
  if (is.null(rownames(mu)) || is.null(rownames(var))) {
    stop("set embedding must have rownames")
  }
  if (nrow(mu) != nrow(var)) stop("set_mu and set_var must have the same number of rows")

  genes <- intersect(names(gene_stats), rownames(gene_emb))
  if (length(genes) < 2) {
    stop("not enough genes overlap between gene_stats and gene_embedding")
  }
  gene_emb <- gene_emb[genes, , drop = FALSE]
  # 与原 api.R 一致：按 embedding 行序重排，保证 ES 与置换同序
  gene_stats <- gene_stats[rownames(gene_emb)]
  storage.mode(gene_stats) <- "double"
  genes <- rownames(gene_emb)

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

  # 先子集 mu/var，再用 1..n_local 分批，避免原脚本在 sets!=NULL 且触发分批时
  # 用「原始行号」去索引已子集矩阵的潜在问题（合回时建议一并修正）。
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
    .proto_weighted_test_sets(
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
      perm_mat = perm_mat
    )
  })

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out$p.adjust <- stats::p.adjust(out$pvalue, method = "BH")
  out
}

# ---- Weighted-ORA：二值 y_g --------------------------------------------------

#' Weighted-ORA（原型）：soft membership 加权超几何风格检验
#'
#' 输入差异基因二值标签 y_g ∈ {0,1}，计算 T_S = sum_g w_gS * y_g，
#' 对 y 做置换得到经验 p，全集合并后 BH-FDR。
#'
#' @param y 命名数值向量，取值只能是 0/1（允许数值型）。
#' @param x gsemb_embedding 对象。
#' @param sets 可选，只测部分 set ID。
#' @param gene_sets 命名 list；报 setSize；restrict_to_members=TRUE 时必填。
#' @param restrict_to_members 是否仅成员参与 soft membership（默认 FALSE）。
#' @param ... 其余传给驱动：score, temperature, nperm, alternative, seed, eps, top_genes。
#' @return data.frame：ID, ES, z, pvalue, p.adjust, setSize, core_enrichment。
gsemb_weighted_ora <- function(y,
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
                               top_genes = 30) {
  if (!is.numeric(y) || is.null(names(y))) {
    stop("y must be a named numeric vector (0/1 labels)")
  }
  # 允许浮点 0/1，拒绝其它值；保留 names（as.numeric 会丢掉名字）
  y_names <- names(y)
  y_num <- as.numeric(y)
  names(y_num) <- y_names
  uniq <- unique(y_num[is.finite(y_num)])
  if (length(uniq) == 0L || !all(uniq %in% c(0, 1))) {
    stop("Weighted-ORA 要求 y 仅含 0 和 1（差异基因二值标签）")
  }
  storage.mode(y_num) <- "double"

  .proto_weighted_enrichment_driver(
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
    top_genes = top_genes
  )
}

# ---- Weighted-GSEA：连续 r_g -------------------------------------------------

#' Weighted-GSEA（原型）：soft membership 加权连续富集
#'
#' 输入连续排序统计量 r_g（如 t、logFC），计算 T_S = sum_g w_gS * r_g。
#' restrict_to_members=FALSE 时应与包内 gsemb_embedding_enrichment 对齐。
#'
#' @param r 命名数值向量（连续基因统计量）。
#' @param x gsemb_embedding 对象。
#' @param sets 可选 set ID。
#' @param gene_sets 命名 list。
#' @param restrict_to_members 是否仅成员参与 soft membership（默认 FALSE）。
#' @param ... 其余同 Weighted-ORA。
#' @return data.frame：ID, ES, z, pvalue, p.adjust, setSize, core_enrichment。
gsemb_weighted_gsea <- function(r,
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
                                top_genes = 30) {
  if (!is.numeric(r) || is.null(names(r))) {
    stop("r must be a named numeric vector (continuous gene-level statistics)")
  }

  .proto_weighted_enrichment_driver(
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
    top_genes = top_genes
  )
}

# ---- 玩具演示（默认不执行）--------------------------------------------------

if (FALSE) {
  # 在包根目录执行：
  #   devtools::load_all(".")
  #   source("scripts/weighted_ora_gsea_prototype.R", encoding = "UTF-8")

  set.seed(1)
  # 列名需为 node1/node2（与 gsemb_fit / gsemb_build_graph 默认一致）
  edges <- data.frame(
    node1 = c("A", "A", "B", "C", "D"),
    node2 = c("B", "C", "C", "D", "E"),
    weight = c(1, 1, 1, 1, 1)
  )
  gene_sets <- list(S1 = c("A", "B", "C"), S2 = c("D", "E"))
  fit <- gsemb_fit(
    edges, gene_sets,
    weight = "weight", method = "svd",
    dim = 2, k = 2, max_iter = 30, epochs = 3
  )

  # 连续统计量 → Weighted-GSEA
  r <- c(A = 2, B = 1, C = 0.5, D = -1, E = -2)
  gsea_all <- gsemb_weighted_gsea(
    r, fit,
    gene_sets = gene_sets,
    restrict_to_members = FALSE,
    nperm = 50, seed = 1, top_genes = 3,
    alternative = "greater"
  )
  gsea_mem <- gsemb_weighted_gsea(
    r, fit,
    gene_sets = gene_sets,
    restrict_to_members = TRUE,
    nperm = 50, seed = 1, top_genes = 3,
    alternative = "greater"
  )

  # 与原函数对照（restrict_to_members=FALSE）
  legacy <- gsemb_embedding_enrichment(
    gene_stats = r, x = fit, gene_sets = gene_sets,
    nperm = 50, seed = 1, top_genes = 3, alternative = "greater"
  )
  cat("GSEA vs legacy max |ES| diff:",
      max(abs(gsea_all$ES - legacy$ES[match(gsea_all$ID, legacy$ID)])), "\n")

  # 二值标签 → Weighted-ORA
  y <- c(A = 1, B = 1, C = 0, D = 0, E = 0)
  ora <- gsemb_weighted_ora(
    y, fit,
    gene_sets = gene_sets,
    restrict_to_members = TRUE,
    nperm = 50, seed = 1, top_genes = 3,
    alternative = "greater"
  )

  print(gsea_all)
  print(gsea_mem)
  print(ora)

  # 可选：真实 fit
  # load("E:/genesetembedding/step3/fit.RData")
  # load("E:/genesetembedding/step3/reactome_gene_sets.RData")
  # ...
}

message(
  "已加载 Weighted-ORA/GSEA 原型：gsemb_weighted_ora() / gsemb_weighted_gsea()；",
  "原 R/api.R 未改动。详见脚本文首合回地图。"
)
