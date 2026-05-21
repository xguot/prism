#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Project the non-parametric imputation estimates onto the structural
 * manifold defined by the target covariance while preserving fidelity to
 * the original imputed values. The augmented loss is:
 *
 *   L(X) = (1/2) ||X - X_imp||_F^2  +  (lambda/2) ||cov(X) - Sigma_target||_F^2
 *
 * The gradient used at each step is:
 *
 *   grad = (X - X_imp)  +  (2 * lambda / (n - 1)) * X_tilde * (cov(X) - Sigma_target)
 *
 * where X_tilde is the column-centered data matrix. Centering is required
 * because the covariance gradient is defined with respect to mean-centred
 * variables; using raw (uncentered) values would leak location shifts into
 * the covariance projection.
 *
 * Parameters:
 *   X_imp        n x p initial imputation matrix (anchor for fidelity term)
 *   Sigma_target p x p target covariance matrix (must be positive semidefinite)
 *   lambda       trade-off weight on the covariance constraint
 *   lr           learning rate for gradient descent
 *   max_iter     maximum number of iterations
 */
// [[Rcpp::export]]
arma::mat constrain_covariance(arma::mat X_imp, arma::mat Sigma_target,
                               double lambda, double lr, int max_iter)
{
  int n = X_imp.n_rows;
  arma::mat X_opt  = X_imp;
  arma::mat X_orig = X_imp;   // anchor for fidelity penalty
  arma::mat X_centered, Sigma_curr, grad_fidelity, grad_cov, Sigma_new;

  for (int i = 0; i < max_iter; i++) {

    // ---- 1.  centre the data (required for correct covariance gradient) ----
    X_centered = X_opt.each_row() - arma::mean(X_opt, 0);

    // ---- 2.  current covariance estimate ----
    Sigma_curr = arma::cov(X_opt);

    // ---- 3.  fidelity gradient: keeps X close to the original imputation ----
    grad_fidelity = X_opt - X_orig;

    // ---- 4.  covariance-constraint gradient ----
    // d/dX ||cov(X) - Sigma||^2 = (4/(n-1)) * X_tilde * (cov(X) - Sigma)
    // With the (lambda/2) weight from the loss: multiply by (2*lambda/(n-1))
    grad_cov =
        (2.0 * lambda / (n - 1)) * X_centered * (Sigma_curr - Sigma_target);

    // ---- 5.  gradient descent step ----
    X_opt = X_opt - lr * (grad_fidelity + grad_cov);

    // ---- 6.  guard against numerical blow-up ----
    if (X_opt.has_nan() || X_opt.has_inf()) {
      Rcpp::stop(
          "Numerical instability in covariance projection. "
          "The target covariance matrix may be ill-conditioned or not "
          "positive semidefinite. Try reducing the learning rate, or verify "
          "the target matrix is valid before calling constrain_covariance().");
    }

    // ---- 7.  convergence check (post-update) ----
    Sigma_new = arma::cov(X_opt);
    if (arma::norm(Sigma_new - Sigma_target, "fro") < 1e-6)
      break;
  }

  return X_opt;
}
