helper_path <- testthat::test_path("..", "..", "scripts", "real_cohort", "02_sle_landmark_helpers.R")
source(helper_path, local = TRUE)

test_that("map_ensp_to_symbols remaps STRING IDs and rejects missing IDs", {
  info <- data.frame(
    string_protein_id = c("9606.ENSP1", "9606.ENSP2"),
    preferred_name = c("AAA", "BBB"),
    stringsAsFactors = FALSE
  )
  expect_equal(map_ensp_to_symbols(c("9606.ENSP2", "9606.ENSP1"), info), c("BBB", "AAA"))
  expect_error(map_ensp_to_symbols("9606.MISSING", info), "Unmapped")
})

test_that("remap_adj_to_symbols renames nodes for DEG/Reactome overlap", {
  adj <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(2, 1),
    x = 1,
    dims = c(2, 2),
    dimnames = list(c("9606.ENSP1", "9606.ENSP2"), c("9606.ENSP1", "9606.ENSP2"))
  )
  info <- data.frame(
    string_protein_id = c("9606.ENSP1", "9606.ENSP2"),
    preferred_name = c("AAA", "BBB"),
    stringsAsFactors = FALSE
  )
  out <- remap_adj_to_symbols(adj, info)
  expect_equal(rownames(out), c("AAA", "BBB"))
  expect_equal(colnames(out), c("AAA", "BBB"))
})

test_that("prepare_ranked_stats keeps the strongest finite statistic per gene", {
  deg <- data.frame(
    gene_symbol = c("A", "A", "B", "", "C"),
    log2FoldChange = c(1, -3, 2, 4, NA_real_),
    stringsAsFactors = FALSE
  )

  stats <- prepare_ranked_stats(deg)

  expect_equal(stats, c(A = -3, B = 2))
})

test_that("build_fit_from_landmarks uses exactly the supplied landmarks", {
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 1),
    j = c(2, 1, 3, 2, 4, 3, 1, 4),
    x = 1,
    dims = c(4, 4),
    dimnames = list(LETTERS[1:4], LETTERS[1:4])
  )
  sets <- list(S1 = c("A", "B"), S2 = c("C", "D"))

  fit <- build_fit_from_landmarks(
    adj = adj,
    gene_sets = sets,
    landmarks = c("A", "C"),
    dim = 1L,
    alpha = 0.5,
    seed = 7L
  )

  expect_s3_class(fit, "gsemb_embedding")
  expect_identical(fit$landmarks, c("A", "C"))
  expect_equal(rownames(fit$gene_embedding), LETTERS[1:4])
  expect_equal(rownames(fit$set_mu), c("S1", "S2"))
})

test_that("as_comparable_enrichment keeps only testable rows and standard columns", {
  x <- data.frame(
    ID = c("P1", "P2"),
    ES = c(2, 1),
    pvalue = c(0.01, NA_real_),
    p.adjust = c(0.02, NA_real_),
    status = c("ok", "untestable_no_overlap"),
    stringsAsFactors = FALSE
  )

  out <- as_comparable_enrichment(x)

  expect_equal(out$ID, "P1")
  expect_true(all(c("ID", "pvalue", "p.adjust", "ES", "status") %in% names(out)))
})

test_that("summarize_method_reproducibility averages pairwise metrics by method", {
  x <- data.frame(
    version = rep(c("weighted_degree_high", "weighted_kmedoids_high"), each = 2),
    pathway_mode = "top20",
    jaccard = c(0.2, 0.4, 0.5, 0.7),
    overlap_coefficient = c(0.3, 0.5, 0.6, 0.8),
    rank_spearman = c(0.1, 0.3, 0.4, 0.6),
    stringsAsFactors = FALSE
  )

  out <- summarize_method_reproducibility(x, route = "weighted_gsea")

  km <- out[out$method == "weighted_kmedoids_high", ]
  expect_equal(km$mean_jaccard, 0.6)
  expect_equal(km$mean_rank_spearman, 0.5)
  expect_identical(km$route, "weighted_gsea")
})
