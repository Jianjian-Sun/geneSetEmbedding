# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`geneSetEmbedding` is an R package that transforms PPI networks into learnable gene embeddings via random-walk diffusion, then represents gene sets as diagonal Gaussians in the embedding space. The package supports both fast SVD-based and GPU-accelerated torch-based embedding methods.

## Build & Test Commands

```bash
# Build package
R CMD build .

# Check package (includes examples and tests)
R CMD check --no-manual geneSetEmbedding_*.tar.gz

# Run tests only
Rscript -e "testthat::test_local('tests/testthat/')"

# Run roxygen2 documentation rebuild
Rscript -e "roxygen2::roxygenise()"

# Quick load for interactive development
Rscript -e "devtools::load_all('.')"
```

## Architecture

### Core Data Flow

```
PPI Edge List → adjacency matrix (dgCMatrix)
               → transition matrix (column/row stochastic)
               → landmark selection (degree-based or random)
               → RWR diffusion (n_nodes × n_landmarks feature matrix)
               → gene embedding (SVD / torch autoencoder / Set2Gaussian torch)
               → set Gaussians (μ, σ² per set from member gene embeddings)
               → gene similarity / set similarity / enrichment analysis
```

### S3 Class: `gsemb_embedding`

Returned by `gsemb_fit()` with fields:
- `adj`: sparse adjacency matrix (dgCMatrix)
- `gene_embedding`: n_genes × dim numeric matrix
- `set_mu`, `set_var`: n_sets × dim numeric matrices
- `landmarks`: character vector of landmark node IDs
- `method`: "svd", "torch_autoencoder", or "set2gaussian_torch"
- `losses`: training loss history (torch methods only)

### Key Files

- **R/api.R**: High-level functions — `gsemb_fit`, `gsemb_gene_similarity`, `gsemb_set_similarity`, `gsemb_embedding_enrichment`, `gsemb_concise_gene_sets`, `gsemb_calculate_all_similarities`, `gsemb_cluster_similarity`
- **R/diffusion.R**: RWR, multi-seed diffusion (`gsemb_rwr`, `gsemb_diffuse_seeds`, `gsemb_compute_node_landmark_features`, `gsemb_compute_set_landmark_features`)
- **R/embedding.R**: Gene embedding via SVD or torch autoencoder; Gaussian fitting from gene set members
- **R/similarity.R**: Cosine similarity for genes; diagonal Gaussian distances (W2, symmetric KL) for sets
- **R/concise.R**: Gene-to-set scoring and concise gene set selection
- **R/set2gaussian_torch.R**: End-to-end Set2Gaussian training with torch
- **R/graph.R**: Graph building from edge lists, transition matrix, landmark selection
- **R/compat_enrichit.R**: Compatibility wrappers (aliases) for enrichit package integration
- **R/00-utils.R**: Shared utilities — `as_numeric_matrix`, `validate_gene_sets`, `require_torch`

### Three Embedding Methods

1. **SVD** (default, `method="svd"`) — direct truncated SVD on diffusion features; fast, CPU-only
2. **torch_autoencoder** (`method="torch_autoencoder"`) — linear encoder/decoder trained with MSE; GPU-accelerated
3. **set2gaussian_torch** (`method="set2gaussian_torch"`) — end-to-end model matching softmax distributions over landmark nodes; most expressive, GPU-accelerated

## Known Issues & Optimization Opportunities

### Fixed Issues

- **diffusion.R:96** — `gsemb_diffuse_seeds` convergence check used `tol * ncol(P)` (too loose); fixed to `tol`.
- **test-set2gaussian-torch.R** — Added `skip_if(!cuda_is_available())` guard for Lantern libtorch environment check.
- **similarity.R** — `gsemb_set_gaussian_distance`: both `w2` and `sym_kl` are now fully vectorized (no double R loop). Benchmarks (n_sets=300, d=64): W2 ~610 ms, sym_KL ~2100 ms.
- **concise.R** — `gsemb_gene_to_set_score` fully vectorized using tiled broadcasting. `gsemb_make_concise_gene_sets` pre-computes all scores in one call when `restrict_to_members=TRUE`. Benchmarks (n_genes=5000, n_sets=300, d=64): gene-to-set ~1.6s batch, concise ~145ms.

### Remaining Opportunities

- **api.R:585-605** — `gsemb_embedding_enrichment` permutation loop (single-threaded R, O(nperm × n_sets)). Consider `Rcpp` for batch permutation + matrix multiply, or `future.apply` for parallel permutations.
- **embedding.R:158** — `apply(X, 2, stats::var)` calls R's `var()` per column; `matrixStats::colVars()` (or single-pass Rcpp) would be faster.
- **concise.R `softmax_mass` select** — The `exp()` + `cumsum()` loop per set is fine for small sets but could be vectorized across all sets simultaneously.

### Rcpp Implementation (src/gaussian_distance.cpp)

Provides `w2_distance()` and `sym_kl_distance()` via RcppArmadillo + OpenMP.
- **W2**: uses Arma BLAS `repmat` + matrix products — same speed as pure-R (BLAS single-threaded in this environment).
- **sym_KL**: BLAS for ratio matrices + OpenMP parallel loop for mahal terms. Scales to 1000 sets (W2 ~7s, KL ~24s on 96-core).
- Fallback pure-R implementation present if Rcpp unavailable.

### RcppArmadillo Build Requirements

- `DESCRIPTION`: `Imports: Rcpp, RcppArmadillo`
- `NAMESPACE`: `useDynLib(geneSetEmbedding, .registration = TRUE)`
- `src/Makevars`: must include both Rcpp and RcppArmadillo include paths, `-fopenmp -lgomp` linker flags. The Arma BLAS gemm calls only parallelize when R itself is compiled with OpenMP (not the case here).

## Rcpp Integration Notes

When adding Rcpp code:
1. Create `src/` directory and C++ source files
2. Add `LinkingTo: Rcpp` and `Imports: Rcpp` to DESCRIPTION
3. Use `.onLoad` in `R/zzz.R` to call `Rcpp::sourceCpp()` if needed, or use `@useDynLib geneSetEmbedding` in NAMESPACE
4. Matrix operations in C++: use `RcppArmadillo` for linear algebra or `RcppEigen` for sparse matrices
5. The package already uses `Matrix::sparseMatrix` extensively — keep sparse representation through Rcpp interfaces

## Dependencies

- **Hard** (Imports): `Matrix`, `stats`
- **Optional** (Suggests): `torch` (GPU training), `testthat` (tests)
- R >= 4.1.0