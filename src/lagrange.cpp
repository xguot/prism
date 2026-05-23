#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Project the non-parametric imputation estimates onto the structural
 * manifold defined by the target covariance while preserving fidelity to
 * the original imputed values.  Only entries marked in the missingness
 * mask are updated; observed data is frozen.
 *
 * Augmented loss (applied to missing entries only):
 *
 *   L(X) = (1/2) ||M .* (X - X_imp)||_F^2
 *        + (lambda/2) ||cov(X) - Sigma_target||_F^2
 *
 * where M is the 0/1 missingness mask (1 = originally missing) and .* is
 * element-wise multiplication.
 *
 * The gradient used at each step (masked so only missing positions move):
 *
 *   grad = M .* (X - X_imp)  +  M .* [ 2*lambda * X_tilde * (cov(X) - Sigma_target) ]
 *
 * Note: the covariance gradient does NOT divide by (n-1).  This keeps the
 * gradient scale comparable to the fidelity term regardless of sample size,
 * preventing the constraint from asymptotically vanishing for large N.
 * Adjust lambda up or down to control the fidelity/constraint trade-off.
 *
 * Parameters:
 *   X_imp        n x p initial imputation matrix
 *   mask         n x p 0/1 missingness indicator (1 = originally missing)
 *   Sigma_target p x p target covariance (must be positive semidefinite)
 *   lambda       per-observation penalty on covariance deviation
 *   lr           learning rate
 *   max_iter     maximum iterations
 *   tol          convergence tolerance (Frobenius norm)
 */
// [[Rcpp::export]]
arma::mat constrain_covariance(arma::mat X_imp, arma::mat mask,
                               arma::mat Sigma_target, double lambda,
                               double lr, int max_iter, double tol)
{
  int n = X_imp.n_rows;
  arma::mat X_opt  = X_imp;
  arma::mat X_orig = X_imp;
  arma::mat X_centered, Sigma_curr, grad_fidelity, grad_cov, Sigma_new;

  for (int i = 0; i < max_iter; i++) {

    /* centre the data */
    X_centered = X_opt.each_row() - arma::mean(X_opt, 0);

    /* current covariance */
    Sigma_curr = arma::cov(X_opt);

    /* fidelity gradient (only non-zero where originally missing) */
    grad_fidelity = mask % (X_opt - X_orig);

    /* covariance-constraint gradient (un-normalised: no 1/(n-1)) */
    /*
     * The augmented loss includes (lambda/2) * ||cov(X) - Sigma_target||_F^2.
     * The gradient of ||cov(X) - Sigma_target||_F^2 w.r.t. X, without the
     * 1/(n-1) normalisation, is 4 * X_tilde * (cov(X) - Sigma_target).
     * With the (lambda/2) factor applied: 2 * lambda * X_tilde * (cov(X) - Sigma_target).
     */
    grad_cov = 2.0 * lambda * X_centered * (Sigma_curr - Sigma_target);

    /* gradient step (masked: only update originally-missing cells) */
    X_opt = X_opt - lr * (mask % (grad_fidelity + grad_cov));

    /* guard against numerical blow-up */
    if (X_opt.has_nan() || X_opt.has_inf()) {
      Rcpp::stop(
          "Numerical instability in covariance projection. "
          "The target covariance matrix may be ill-conditioned or not "
          "positive semidefinite. Try reducing the learning rate, or verify "
          "the target matrix is valid before calling constrain_covariance().");
    }

    /* convergence check (post-update) */
    Sigma_new = arma::cov(X_opt);
    if (arma::norm(Sigma_new - Sigma_target, "fro") < tol)
      break;
  }

  return X_opt;
}