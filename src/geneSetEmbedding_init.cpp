// geneSetEmbedding_init.cpp
// Custom R_init that registers C++ native functions without RcppExports.
// Do NOT add // [[Rcpp::export]] markers here.

#define ARMA_USE_OPENMP
#include <RcppArmadillo.h>

#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ===========================================================================
// w2_distance — diagonal 2-Wasserstein distance (BLAS)
// ===========================================================================

RcppExport SEXP _geneSetEmbedding_w2_distance(SEXP mu, SEXP var, SEXP mu2, SEXP var2) {
    BEGIN_RCPP
    arma::mat mu_in   = as<arma::mat>(mu);
    arma::mat var_in  = as<arma::mat>(var);
    arma::mat mu2_in  = as<arma::mat>(mu2);
    arma::mat var2_in = as<arma::mat>(var2);

    const uword n1 = mu_in.n_rows;
    const uword n2 = mu2_in.n_rows;

    vec mu_sq  = sum(mu_in  % mu_in,  1);
    vec mu2_sq = sum(mu2_in % mu2_in, 1);

    mat dmu2 = repmat(mu_sq, 1, n2) + repmat(mu2_sq.t(), n1, 1) - 2.0 * mu_in * mu2_in.t();

    mat sv  = sqrt(var_in);
    mat sv2 = sqrt(var2_in);
    mat var_term = repmat(sum(sv  % sv,  1), 1, n2) +
                   repmat(sum(sv2 % sv2, 1).t(), n1, 1) -
                   2.0 * sv * sv2.t();

    return Rcpp::wrap(dmu2 + var_term);
    END_RCPP
}

// ===========================================================================
// sym_kl_distance — symmetric KL divergence between diagonal Gaussians
// Implements: sym_KL(i,j) = 0.5 * (sum_k a[i,k]/b[j,k] + sum_k b[j,k]/a[i,k]
//                                        + mahal_ij[i,j] + mahal_ji[i,j])
// where mahal_ij[i,j] = sum_k (mu[i,k]-mu2[j,k])^2 / b[j,k]
//       mahal_ji[i,j] = sum_k (mu2[i,k]-mu[j,k])^2 / a[j,k]
// Uses outer products + weighted matrix multiplication (BLAS).
// ===========================================================================

RcppExport SEXP _geneSetEmbedding_sym_kl_distance(SEXP mu, SEXP var, SEXP mu2, SEXP var2) {
    BEGIN_RCPP
    arma::mat mu_in   = as<arma::mat>(mu);
    arma::mat var_in  = as<arma::mat>(var);
    arma::mat mu2_in  = as<arma::mat>(mu2);
    arma::mat var2_in = as<arma::mat>(var2);

    const uword n1 = mu_in.n_rows;
    const uword n2 = mu2_in.n_rows;
    const uword d  = mu_in.n_cols;

    // Ratio terms: ratio_mat[i,j] = sum_k var[i,k] * inv_var2[j,k]
    //              ratio_mat2[j,i] = sum_k var2[j,k] * inv_var[i,k]
    mat inv_var  = 1.0 / var_in;
    mat inv_var2 = 1.0 / var2_in;
    mat ratio_mat  = var_in  * inv_var2.t();  // n1 x n2
    mat ratio_mat2 = var2_in * inv_var.t();   // n2 x n1

    // Row sums of mu (n1 x 1) and mu2 (n2 x 1)
    vec mu_sq  = sum(mu_in  % mu_in,  1);
    vec mu2_sq = sum(mu2_in % mu2_in, 1);

    // Row sums of inv_var and inv_var2
    vec iv_rsum  = sum(inv_var,  1);   // n1
    vec iv2_rsum = sum(inv_var2, 1);   // n2

    // mahal_ij[i,j] = mu_sq[i]*iv2_rsum[j] + mu2_sq[j]*iv2_rsum[j]
    //                 - 2 * sum_k mu[i,k] * mu2[j,k] * inv_var2[j,k]
    // = outer(mu_sq, iv2_rsum) + outer(mu2_sq, iv2_rsum) - 2 * S
    // where S[i,j] = sum_k mu[i,k] * mu2[j,k] * inv_var2[j,k]
    //
    // S = mu_in * (mu2_in .* inv_var2).t()   (BLAS matmul)
    mat mu2_w = mu2_in;                         // n2 x d
    for (uword j = 0; j < n2; j++) mu2_w.row(j) %= inv_var2.row(j);  // n2 x d
    mat S = mu_in * mu2_w.t();                  // n1 x n2
    mat mahal_ij = repmat(mu_sq, 1, n2) + repmat(mu2_sq.t(), n1, 1) - 2.0 * S;

    // mahal_ji[i,j] = mu2_sq[i]*iv_rsum[j] + mu_sq[j]*iv_rsum[j]
    //                 - 2 * sum_k mu2[i,k] * mu[j,k] * inv_var[j,k]
    mat mu_w = mu_in;                           // n1 x d
    for (uword i = 0; i < n1; i++) mu_w.row(i) %= inv_var.row(i);   // n1 x d
    mat S2 = mu2_in * mu_w.t();                 // n2 x n1
    mat mahal_ji = repmat(mu2_sq.t(), n1, 1) + repmat(mu_sq, 1, n2) - 2.0 * S2.t();

    return Rcpp::wrap(0.5 * (ratio_mat + ratio_mat2.t() + mahal_ij + mahal_ji));
    END_RCPP
}

// ===========================================================================
// rcpp_enrichment_permutations — parallel Fisher-Yates permutations
// ===========================================================================

RcppExport SEXP _geneSetEmbedding_rcpp_enrichment_permutations(SEXP W, SEXP stats, SEXP n_perm, SEXP seed) {
    BEGIN_RCPP
    arma::mat W_in      = as<arma::mat>(W);
    arma::vec stats_in  = as<arma::vec>(stats);
    int n_perm_in = as<int>(n_perm);
    int seed_in   = as<int>(seed);

    const uword n_gene = W_in.n_rows;
    const uword n_set  = W_in.n_cols;
    mat null_scores(n_perm_in, n_set, fill::zeros);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < n_perm_in; b++) {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
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

// ===========================================================================
// R_init hook — register all native routines
// ===========================================================================

static const R_CallMethodDef CallEntries[] = {
    {"_geneSetEmbedding_w2_distance",                (DL_FUNC) &_geneSetEmbedding_w2_distance,                4},
    {"_geneSetEmbedding_sym_kl_distance",             (DL_FUNC) &_geneSetEmbedding_sym_kl_distance,             4},
    {"_geneSetEmbedding_rcpp_enrichment_permutations", (DL_FUNC) &_geneSetEmbedding_rcpp_enrichment_permutations, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_geneSetEmbedding(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
