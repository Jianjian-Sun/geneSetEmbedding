# BRCA real-cohort benchmark for three Landmark families.
#
# Same design as scripts/real_cohort/01_run_sle_landmark_benchmark.R:
#   fixed Landmarks -> fit -> Weighted-GSEA (primary)
#   fixed Landmarks -> fit -> method-specific Concise sets -> ORA (sensitivity)
#
# Run from R or RStudio:
#   source("scripts/real_cohort/04_run_brca_landmark_benchmark.R")
#
# Outputs go to a new folder under test_results/real_cohort/.
# Fits/Concise sets are disease-independent and reused from the SLE run when present.

# =============================================================================
# 1. Parameters
# =============================================================================
workspace_root <- "E:/genesetembedding"
pkg_root <- file.path(
  workspace_root,
  "geneSetEmbedding/.worktrees/weighted-degree-landmarks"
)
cohorts <- c("GSE109169", "GSE29044", "TCGA_BRCA")

ppi_score_cutoff <- 700
k <- 128L
embedding_dim <- 64L
alpha <- 0.5
tol <- 1e-10
max_iter <- 200L
seed <- 20260919L

gsea_nperm <- 1000L
gsea_alternative <- "greater"
gsea_temperature <- 1
gsea_degree_beta <- 0
gsea_restrict_to_members <- FALSE

concise_select <- "softmax_mass"
concise_mass <- 0.8
concise_temperature <- 1
concise_min_size <- 5L
concise_max_size <- 500L
ora_padj_cutoff <- 0.05
ora_lfc_cutoff <- 1

top_k_values <- c(10L, 20L, 50L)
overwrite <- FALSE

run_label <- "brca_landmark_benchmark_v1"
out_root <- file.path(pkg_root, "test_results", "real_cohort", run_label)
# Same landmarks/Reactome/RWR/SVD as SLE -> reuse cached fits/concise when available.
shared_cache_root <- file.path(
  pkg_root, "test_results", "real_cohort", "sle_landmark_benchmark_v1"
)

# =============================================================================
# 2. Existing inputs
# =============================================================================
high_graph_rds <- file.path(pkg_root, "test_results", "high_conf_graph.rds")
betweenness_rds <- file.path(
  pkg_root, "test_results", "betweenness", "betweenness_landmarks.rds"
)
kmedoids_rds <- file.path(
  pkg_root, "test_results", "kmedoids", "kmedoids_landmarks.rds"
)
string_info_txt <- file.path(
  workspace_root, "step1", "STRING", "9606.protein.info.v12.0.txt"
)
reactome_rdata <- file.path(workspace_root, "step3", "reactome_gene_sets.RData")
reactome_csv <- file.path(workspace_root, "step1", "msigdbr", "reactome_df.csv")
deg_dir <- file.path(workspace_root, "step3", "BRCA", "BRCA_results")
deg_files <- stats::setNames(
  file.path(deg_dir, paste0(cohorts, "_deg_all.csv")),
  cohorts
)
modular_dir <- file.path(workspace_root, "step3", "step3_modular")

required_files <- c(
  high_graph_rds, betweenness_rds, kmedoids_rds, string_info_txt,
  reactome_rdata, reactome_csv, deg_files,
  file.path(modular_dir, "00_helpers.R"),
  file.path(modular_dir, "01_ora_from_deg.R"),
  file.path(modular_dir, "02_compare_ora.R"),
  file.path(pkg_root, "scripts", "real_cohort", "02_sle_landmark_helpers.R")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required files are missing:\n  ", paste(missing_files, collapse = "\n  "), call. = FALSE)
}

# =============================================================================
# 3. Load package + helpers
# =============================================================================
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to load this worktree.", call. = FALSE)
}
devtools::load_all(pkg_root, quiet = TRUE)
source(file.path(pkg_root, "scripts", "real_cohort", "02_sle_landmark_helpers.R"))
.step3_modular_dir <- modular_dir
source(file.path(modular_dir, "00_helpers.R"))
source(file.path(modular_dir, "01_ora_from_deg.R"))
source(file.path(modular_dir, "02_compare_ora.R"))

load(reactome_rdata)
if (!exists("reactome_gene_sets", inherits = FALSE)) {
  stop("reactome_gene_sets.RData does not contain 'reactome_gene_sets'.", call. = FALSE)
}
reactome_df <- utils::read.csv(reactome_csv, stringsAsFactors = FALSE)
protein_info <- load_string_symbol_map(string_info_txt)
adj_high <- remap_adj_to_symbols(readRDS(high_graph_rds), protein_info)

betweenness_obj <- readRDS(betweenness_rds)
kmedoids_obj <- readRDS(kmedoids_rds)
landmarks <- list(
  weighted_degree_high = gsemb_select_landmarks(
    adj_high, k = k, method = "weighted_degree", seed = seed
  ),
  weighted_betweenness_high = map_ensp_to_symbols(
    betweenness_obj$landmarks$weighted_betweenness_high, protein_info
  ),
  weighted_kmedoids_high = map_ensp_to_symbols(
    kmedoids_obj$landmarks$weighted_kmedoids_high, protein_info
  )
)
bad_sizes <- names(landmarks)[vapply(landmarks, length, integer(1)) != k]
if (length(bad_sizes)) {
  stop("Landmark sets do not contain k nodes: ", paste(bad_sizes, collapse = ", "), call. = FALSE)
}

# =============================================================================
# 4. Build/reuse one fit per Landmark method
# =============================================================================
fits_dir <- file.path(out_root, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
shared_fits_dir <- file.path(shared_cache_root, "fits")
fits <- list()
for (method_name in names(landmarks)) {
  fit_path <- file.path(fits_dir, paste0("fit_", method_name, ".rds"))
  shared_fit_path <- file.path(shared_fits_dir, paste0("fit_", method_name, ".rds"))
  if (file.exists(fit_path) && !overwrite) {
    message("Reusing fit: ", fit_path)
    fits[[method_name]] <- readRDS(fit_path)
  } else if (file.exists(shared_fit_path) && !overwrite) {
    message("Reusing shared fit: ", shared_fit_path)
    fits[[method_name]] <- readRDS(shared_fit_path)
    saveRDS(fits[[method_name]], fit_path)
  } else {
    message("Building fit: ", method_name)
    fits[[method_name]] <- build_fit_from_landmarks(
      adj = adj_high,
      gene_sets = reactome_gene_sets,
      landmarks = landmarks[[method_name]],
      dim = embedding_dim,
      alpha = alpha,
      tol = tol,
      max_iter = max_iter,
      normalize = "col",
      seed = seed
    )
    saveRDS(fits[[method_name]], fit_path)
  }
}

# =============================================================================
# 5. Route A (primary): Weighted-GSEA
# =============================================================================
gsea_dir <- file.path(out_root, "weighted_gsea")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
gsea_results <- stats::setNames(vector("list", length(fits)), names(fits))

for (method_name in names(fits)) {
  gsea_results[[method_name]] <- stats::setNames(vector("list", length(cohorts)), cohorts)
  for (cohort in cohorts) {
    result_path <- file.path(gsea_dir, paste0(cohort, "_", method_name, "_weighted_gsea.csv"))
    if (file.exists(result_path) && !overwrite) {
      result <- utils::read.csv(result_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      deg <- utils::read.csv(deg_files[[cohort]], stringsAsFactors = FALSE, check.names = FALSE)
      ranked_stats <- prepare_ranked_stats(deg)
      overlap <- intersect(names(ranked_stats), rownames(fits[[method_name]]$gene_embedding))
      if (length(overlap) < 100L) {
        stop(cohort, " / ", method_name, " has fewer than 100 genes overlapping the fit.", call. = FALSE)
      }
      message("Weighted-GSEA: ", method_name, " / ", cohort)
      result <- gsemb_weighted_gsea(
        r = ranked_stats,
        x = fits[[method_name]],
        gene_sets = reactome_gene_sets,
        temperature = gsea_temperature,
        nperm = gsea_nperm,
        alternative = gsea_alternative,
        seed = seed,
        restrict_to_members = gsea_restrict_to_members,
        degree_beta = gsea_degree_beta
      )
      result <- as_comparable_enrichment(result)
      utils::write.csv(result, result_path, row.names = FALSE)
    }
    gsea_results[[method_name]][[cohort]] <- result
  }
}

gsea_topk <- compare_ora_at_topk_list(
  ora_files = gsea_results,
  top_k_list = top_k_values,
  padj_cutoff = ora_padj_cutoff,
  cohort_order = cohorts,
  out_dir = file.path(gsea_dir, "cross_cohort"),
  out_prefix = "BRCA_landmark_weighted_gsea"
)
gsea_all_sig <- compare_ora_results(
  ora_files = gsea_results,
  padj_cutoff = ora_padj_cutoff,
  top_k = NULL,
  cohort_order = cohorts,
  out_dir = file.path(gsea_dir, "cross_cohort"),
  out_prefix = "BRCA_landmark_weighted_gsea"
)
gsea_pairs <- rbind(gsea_topk$pairs, gsea_all_sig$pairs)
gsea_summary <- summarize_method_reproducibility(gsea_pairs, route = "weighted_gsea")

# =============================================================================
# 6. Route B (sensitivity): Concise sets -> ORA
# =============================================================================
concise_dir <- file.path(out_root, "concise_ora")
dir.create(concise_dir, recursive = TRUE, showWarnings = FALSE)
shared_concise_dir <- file.path(shared_cache_root, "concise_ora")
concise_sets <- list()
for (method_name in names(fits)) {
  concise_path <- file.path(concise_dir, paste0("concise_sets_", method_name, ".rds"))
  shared_concise_path <- file.path(shared_concise_dir, paste0("concise_sets_", method_name, ".rds"))
  if (file.exists(concise_path) && !overwrite) {
    concise_sets[[method_name]] <- readRDS(concise_path)
  } else if (file.exists(shared_concise_path) && !overwrite) {
    message("Reusing shared Concise sets: ", shared_concise_path)
    concise_sets[[method_name]] <- readRDS(shared_concise_path)
    saveRDS(concise_sets[[method_name]], concise_path)
  } else {
    message("Building Concise sets: ", method_name)
    concise_sets[[method_name]] <- gsemb_concise_gene_sets(
      fits[[method_name]],
      gene_sets = reactome_gene_sets,
      select = concise_select,
      mass = concise_mass,
      temperature = concise_temperature,
      restrict_to_members = TRUE,
      min_size = concise_min_size,
      max_size = concise_max_size
    )
    saveRDS(concise_sets[[method_name]], concise_path)
  }
}

ora_results <- run_ora_from_deg_files(
  deg_files = deg_files,
  gene_sets_versions = concise_sets,
  reactome_df = reactome_df,
  out_dir = file.path(concise_dir, "cohort_results"),
  padj_cutoff = ora_padj_cutoff,
  lfc_cutoff = ora_lfc_cutoff,
  plot_top_n = 20L,
  save = TRUE
)
ora_files <- lapply(names(concise_sets), function(method_name) {
  stats::setNames(lapply(cohorts, function(cohort) {
    ora_results[[cohort]]$ora[[method_name]]
  }), cohorts)
})
names(ora_files) <- names(concise_sets)

ora_topk <- compare_ora_at_topk_list(
  ora_files = ora_files,
  top_k_list = top_k_values,
  padj_cutoff = ora_padj_cutoff,
  cohort_order = cohorts,
  out_dir = file.path(concise_dir, "cross_cohort"),
  out_prefix = "BRCA_landmark_concise_ora"
)
ora_all_sig <- compare_ora_results(
  ora_files = ora_files,
  padj_cutoff = ora_padj_cutoff,
  top_k = NULL,
  cohort_order = cohorts,
  out_dir = file.path(concise_dir, "cross_cohort"),
  out_prefix = "BRCA_landmark_concise_ora"
)
ora_pairs <- rbind(ora_topk$pairs, ora_all_sig$pairs)
ora_summary <- summarize_method_reproducibility(ora_pairs, route = "concise_ora")

# =============================================================================
# 7. Summary + provenance
# =============================================================================
method_summary <- rbind(gsea_summary, ora_summary)
method_summary <- method_summary[order(
  method_summary$route,
  method_summary$pathway_mode,
  -method_summary$mean_rank_spearman,
  -method_summary$mean_jaccard
), , drop = FALSE]
utils::write.csv(method_summary, file.path(out_root, "method_summary.csv"), row.names = FALSE)

write_run_metadata(file.path(out_root, "run_meta.txt"), list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  disease = "BRCA",
  pkg_root = pkg_root,
  cohorts = cohorts,
  deg_dir = deg_dir,
  methods = names(landmarks),
  shared_cache_root = shared_cache_root,
  string_info = string_info_txt,
  ppi_score_cutoff = ppi_score_cutoff,
  k = k,
  embedding_dim = embedding_dim,
  alpha = alpha,
  tol = tol,
  max_iter = max_iter,
  seed = seed,
  gsea_nperm = gsea_nperm,
  gsea_alternative = gsea_alternative,
  gsea_degree_beta = gsea_degree_beta,
  concise_select = concise_select,
  concise_mass = concise_mass,
  ora_padj_cutoff = ora_padj_cutoff,
  ora_lfc_cutoff = ora_lfc_cutoff,
  top_k_values = top_k_values
))

message("Finished. Main summary: ", file.path(out_root, "method_summary.csv"))
print(method_summary)
invisible(list(
  landmarks = landmarks,
  fits = fits,
  weighted_gsea = list(results = gsea_results, pairs = gsea_pairs, summary = gsea_summary),
  concise_ora = list(results = ora_results, pairs = ora_pairs, summary = ora_summary),
  method_summary = method_summary
))
