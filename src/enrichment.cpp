// enrichment.cpp
// Parallel permutation enrichment test using RcppArmadillo + OpenMP.
// Note: no [[Rcpp::export]] — compiled into shared library, called via native symbol.

#define ARMA_USE_OPENMP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

using namespace Rcpp;
using namespace arma;

//' Parallel permutations for gene-set enrichment
//'
//' Runs n_perm Fisher-Yates shuffles of gene_stats and computes
//' W^T %*% perm_stats for each (equivalent to R's crossprod(W, perm_stats)).
//' Outer loop over permutations is parallelized via OpenMP.
//'
//' @param W       n_genes x n_sets weight matrix (column-stochastic)
//' @param stats   n_genes gene-level statistics
//' @param n_perm  number of permutations
//' @param seed    random seed
//' @return        n_perm x n_sets matrix of null enrichment scores
// [[Rcpp::export]]
RcppExport SEXP _geneSetEmbedding_rcpp_enrichment_permutations(SEXP W, SEXP stats, SEXP n_perm, SEXP seed) {
    BEGIN_RCPP
    arma::mat W_in    = as<arma::mat>(W);
    arma::vec stats_in = as<arma::vec>(stats);
    int n_perm_in = as<int>(n_perm);
    int seed_in   = as<int>(seed);

    const uword n_gene = W_in.n_rows;
    const uword n_set  = W_in.n_cols;

    mat null_scores(n_perm_in, n_set, fill::zeros);

    #pragma omp for schedule(static)
    for (int b = 0; b < n_perm_in; b++) {
        int t = omp_get_thread_num();
        uint64_t rng_seed = static_cast<uint64_t>(seed_in) * 1315423911ULL +
                             static_cast<uint64_t>(b) * 2654435761ULL +
                             static_cast<uint64_t>(t) * 671936453ULL;
        std::mt19937 rng(rng_seed);
        std::uniform_int_distribution<uword> draw(0, n_gene - 1);

        uvec idx(n_gene);
        for (uword k = 0; k < n_gene; k++) idx[k] = k;
        for (uword k = n_gene; k > 1; k--) {
            uword r = draw(rng) % k;
            uword tmp = idx[k - 1]; idx[k - 1] = idx[r]; idx[r] = tmp;
        }

        vec perm_stats = stats_in.elem(idx);
        vec es_null = W_in.t() * perm_stats;

        double* out_row = null_scores.memptr() + b * n_set;
        const double* in_row = es_null.memptr();
        for (uword s = 0; s < n_set; s++) out_row[s] = in_row[s];
    }

    return Rcpp::wrap(null_scores);
    END_RCPP
}