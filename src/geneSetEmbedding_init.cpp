// geneSetEmbedding_init.cpp
// Custom R_init that registers C++ native functions without RcppExports.
// Compiled from: gaussian.cpp (w2_distance, sym_kl_distance) and enrichment.cpp (rcpp_enrichment_permutations).

#define ARMA_USE_OPENMP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

using namespace Rcpp;
using namespace arma;

// ===========================================================================
// gaussian.cpp — w2_distance
// ===========================================================================

//' @export
// [[Rcpp::export]]
RcppExport SEXP _geneSetEmbedding_w2_distance(SEXP mu, SEXP var, SEXP mu2, SEXP var2) {
    BEGIN_RCPP
    arma::mat mu_in  = as<arma::mat>(mu);
    arma::mat var_in = as<arma::mat>(var);
    arma::mat mu2_in  = as<arma::mat>(mu2);
    arma::mat var2_in = as<arma::mat>(var2);

    const uword n1 = mu_in.n_rows;
    const uword n2 = mu2_in.n_rows;
    const uword d  = mu_in.n_cols;

    vec mu_sq  = sum(mu_in  % mu_in,  1);
    vec mu2_sq = sum(mu2_in % mu2_in, 1);

    // ||mu[i] - mu2[j]||^2 = repmat(mu_sq,1,n2) + repmat(mu2_sq.t(),n1,1) - 2*mu*mu2.t()
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
// gaussian.cpp — sym_kl_distance
// ===========================================================================

//' @export
// [[Rcpp::export]]
RcppExport SEXP _geneSetEmbedding_sym_kl_distance(SEXP mu, SEXP var, SEXP mu2, SEXP var2) {
    BEGIN_RCPP
    arma::mat mu_in  = as<arma::mat>(mu);
    arma::mat var_in = as<arma::mat>(var);
    arma::mat mu2_in  = as<arma::mat>(mu2);
    arma::mat var2_in = as<arma::mat>(var2);

    const uword n1 = mu_in.n_rows;
    const uword n2 = mu2_in.n_rows;
    const uword d  = mu_in.n_cols;

    mat inv_var  = 1.0 / var_in;
    mat inv_var2 = 1.0 / var2_in;

    mat ratio_mat  = var_in  * inv_var2.t();  // n1 x n2
    mat ratio_mat2 = var2_in * inv_var.t();   // n2 x n1

    vec mu_sq    = sum(mu_in  % mu_in,  1);
    vec mu2_sq   = sum(mu2_in % mu2_in, 1);
    vec iv_rsum  = sum(inv_var,  1);   // n1
    vec iv2_rsum = sum(inv_var2, 1);   // n2

    // mahal_ij = repmat(mu_sq,1,n2)%repmat(iv2_rsum.t(),n1,1) + repmat(mu2_sq.t(),n1,1)%repmat(iv2_rsum.t(),n1,1) - 2*mu%*%(mu2%*%diag(iv2_rsum))
    mat mu2_w = mu2_in.each_col() % inv_var2;   // n2 x d
    mat mahal_ij = repmat(mu_sq,   1, n2) % repmat(iv2_rsum.t(), n1, 1) +
                   repmat(mu2_sq.t(), n1, 1) % repmat(iv2_rsum.t(), n1, 1) -
                   2.0 * mu_in * mu2_w.t();

    mat mu_w = mu_in.each_col() % inv_var;
    mat mahal_ji = repmat(mu2_sq.t(), n1, 1) % repmat(iv_rsum.t(), n1, 1) +
                   repmat(mu_sq,   1, n2) % repmat(iv_rsum.t(), n1, 1) -
                   2.0 * mu2_in * mu_w.t();

    return Rcpp::wrap(0.5 * (ratio_mat + ratio_mat2.t() + mahal_ij + mahal_ji));
    END_RCPP
}

// ===========================================================================
// enrichment.cpp — rcpp_enrichment_permutations (parallel)
// ===========================================================================

//' @export
// [[Rcpp::export]]
RcppExport SEXP _geneSetEmbedding_rcpp_enrichment_permutations(SEXP W, SEXP stats, SEXP n_perm, SEXP seed) {
    BEGIN_RCPP
    arma::mat W_in      = as<arma::mat>(W);
    arma::vec stats_in  = as<arma::vec>(stats);
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

// ===========================================================================
// R_init hook — register all native functions
// ===========================================================================

static const R_CallMethodDef CallEntries[] = {
    {"_geneSetEmbedding_w2_distance",               (DL_FUNC) &_geneSetEmbedding_w2_distance,               4},
    {"_geneSetEmbedding_sym_kl_distance",            (DL_FUNC) &_geneSetEmbedding_sym_kl_distance,            4},
    {"_geneSetEmbedding_rcpp_enrichment_permutations", (DL_FUNC) &_geneSetEmbedding_rcpp_enrichment_permutations, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_geneSetEmbedding(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}