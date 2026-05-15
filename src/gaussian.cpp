// gaussian_distance.cpp
// RcppArmadillo + OpenMP implementation of pairwise Gaussian distance functions.
// Mahal term uses OpenMP parallel loops; all BLAS calls remain available.

#define ARMA_USE_OPENMP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>

using namespace Rcpp;
using namespace arma;

//' Diagonal 2-Wasserstein distance (fully vectorized via BLAS)
//'
//' out[i,j] = ||mu[i] - mu2[j]||^2 + ||sqrt(var[i]) - sqrt(var2[j])||^2
//'
//' @param mu   n1 x d matrix of means
//' @param var  n1 x d matrix of variances (must be > 0)
//' @param mu2  n2 x d matrix of means
//' @param var2 n2 x d matrix of variances (must be > 0)
//' @return     n1 x n2 matrix of squared distances
// [[Rcpp::export]]
arma::mat w2_distance(const arma::mat& mu, const arma::mat& var,
                       const arma::mat& mu2, const arma::mat& var2) {
    const uword n1 = mu.n_rows;
    const uword n2 = mu2.n_rows;
    const uword d  = mu.n_cols;

    vec mu_sq  = sum(mu  % mu,  1);
    vec mu2_sq = sum(mu2 % mu2, 1);

    // ||mu[i] - mu2[j]||^2 = ||mu[i]||^2 + ||mu2[j]||^2 - 2 * mu[i] · mu2[j]
    mat dmu2 = repmat(mu_sq,  1, n2) + repmat(mu2_sq.t(), n1, 1) - 2.0 * mu * mu2.t();

    mat sv  = sqrt(var);
    mat sv2 = sqrt(var2);
    mat var_term = repmat(sum(sv  % sv,  1), 1, n2) +
                   repmat(sum(sv2 % sv2, 1).t(), n1, 1) -
                   2.0 * sv * sv2.t();

    return dmu2 + var_term;
}

//' Symmetric KL divergence between diagonal Gaussians
//'
//' sym_KL(i,j) = 0.5 * ( sum_d var[i,d]/var2[j,d] + sum_d var2[j,d]/var[i,d]
//'                       + mahal_ij[i,j] + mahal_ji[i,j] )
//'
//' mahal_ij[i,j] = sum_d (mu[i,d]-mu2[j,d])^2 / var2[j,d]
//'
//' @param mu   n1 x d matrix of means
//' @param var  n1 x d matrix of variances
//' @param mu2  n2 x d matrix of means
//' @param var2 n2 x d matrix of variances
//' @return     n1 x n2 matrix of symmetric KL divergences
// [[Rcpp::export]]
arma::mat sym_kl_distance(const arma::mat& mu, const arma::mat& var,
                            const arma::mat& mu2, const arma::mat& var2) {
    const uword n1 = mu.n_rows;
    const uword n2 = mu2.n_rows;
    const uword d  = mu.n_cols;

    // Variance ratio matrices via BLAS: O(n1 * n2 * d)
    mat inv_var  = 1.0 / var;
    mat inv_var2 = 1.0 / var2;
    mat ratio_mat  = var  * inv_var2.t();  // n1 x n2
    mat ratio_mat2 = var2 * inv_var.t();   // n2 x n1

    // Pre-compute per-row scalars
    vec mu_sq   = sum(mu  % mu,  1);
    vec mu2_sq  = sum(mu2 % mu2, 1);
    vec iv_rsum  = sum(inv_var,  1);   // n1
    vec iv2_rsum = sum(inv_var2, 1);   // n2

    // Allocate output and mahal buffers
    mat out(n1, n2, fill::zeros);
    mat mahal_ij(n1, n2, fill::zeros);
    mat mahal_ji(n1, n2, fill::zeros);

    // Get raw pointers for OpenMP loop
    const double* mu_ptr   = mu.memptr();
    const double* mu2_ptr  = mu2.memptr();
    const double* iv_ptr   = inv_var.memptr();
    const double* iv2_ptr  = inv_var2.memptr();
    const double* mu_sq_ptr   = mu_sq.memptr();
    const double* mu2_sq_ptr  = mu2_sq.memptr();
    const double* iv_rsum_ptr  = iv_rsum.memptr();
    const double* iv2_rsum_ptr = iv2_rsum.memptr();
    double* m12_ptr = mahal_ij.memptr();
    double* m21_ptr = mahal_ji.memptr();

    // Parallel outer loop over rows of mu (n1) — all inner computations independent
    #pragma omp parallel for schedule(static)
    for (uword i = 0; i < n1; i++) {
        const double* mu_i   = mu_ptr   + i * d;
        const double* iv_i   = iv_ptr   + i * d;
        double mu_sq_i       = mu_sq_ptr[i];
        double iv_rsum_i     = iv_rsum_ptr[i];

        for (uword j = 0; j < n2; j++) {
            const double* mu2_j  = mu2_ptr  + j * d;
            const double* iv2_j  = iv2_ptr  + j * d;
            double mu2_sq_j      = mu2_sq_ptr[j];
            double iv2_rsum_j    = iv2_rsum_ptr[j];

            // mahal_ij[i,j] = sum_d (mu[i,d] - mu2[j,d])^2 / var2[j,d]
            double dot12 = 0.0;
            for (uword k = 0; k < d; k++) {
                double diff = mu_i[k] - mu2_j[k];
                dot12 += diff * diff * iv2_j[k];
            }
            m12_ptr[i + j * n1] = dot12;

            // mahal_ji[i,j] = sum_d (mu2[j,d] - mu[i,d])^2 / var[i,d]
            // Symmetric with i/j swapped: reuse loop over k with mu2[j], mu[i], inv_var[i]
            double dot21 = 0.0;
            for (uword k = 0; k < d; k++) {
                double diff = mu2_j[k] - mu_i[k];
                dot21 += diff * diff * iv_i[k];
            }
            m21_ptr[i + j * n1] = dot21;
        }
    }

    out = 0.5 * (ratio_mat + ratio_mat2.t() + mahal_ij + mahal_ji);
    return out;
}
