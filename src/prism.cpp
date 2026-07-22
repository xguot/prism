#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Project the non-parametric imputation estimates onto the structural
 * manifold defined by the target covariance. Only entries marked in the
 * missingness mask are updated; observed data is frozen.
 *
 * This optimized version focuses on structural covariance preservation
 * using a Lagrangian constraint.
 *
 * Parameters:
 *   X_imp        n x p initial imputation matrix
 *   mask         n x p 0/1 missingness indicator (1 = originally missing)
 *   Sigma_target p x p target covariance (must be positive semidefinite)
 *   lambda       per-observation penalty on covariance deviation
 *   lr           learning rate
 *   max_iter     maximum iterations
 *   tol          convergence tolerance (Frobenius norm)
 *
 * Returns a list with:
 *   X_refined    n x p refined matrix
 *   iterations   number of gradient descent iterations executed
 *   converged    logical; TRUE if Frobenius norm fell below tol
 *   final_frob   final Frobenius distance to target covariance
 */
// [[Rcpp::export]]
Rcpp::List constrain_covariance(const arma::mat& X_imp,
                                const arma::mat& mask,
                                const arma::mat& Sigma_target,
                                double lambda, double lr, int max_iter, double tol) {

  arma::mat X_opt = X_imp;
  int n = X_opt.n_rows;
  arma::mat X_centered, Sigma_curr, grad_cov, update_step;

  /* Strict PSD Guard: Ensure the target is positive semi-definite.
   * Copy to a local mutable matrix — modifying a const reference via
   * const_cast is UB.  Outliers in heavy distributions can occasionally
   * cause precision issues that lead to negative eigenvalues, even if
   * the R-side already projected.
   */
  arma::mat Sigma_fixed = Sigma_target;
  arma::vec eigval_t;
  arma::mat eigvec_t;
  if (arma::eig_sym(eigval_t, eigvec_t, Sigma_fixed)) {
    bool has_neg = false;
    for (uword j = 0; j < eigval_t.n_elem; j++) {
      if (eigval_t(j) < 0) {
        eigval_t(j) = 0;
        has_neg = true;
      }
    }
    if (has_neg) {
      Sigma_fixed = eigvec_t * arma::diagmat(eigval_t) * eigvec_t.t();
    }
  }

  int iter = 0;
  bool converged = false;
  double final_frob = -1.0;

  for (int i = 0; i < max_iter; i++) {
    /* check for user interrupt to allow Esc/Ctrl+C in R */
    Rcpp::checkUserInterrupt();

    /* efficient centering and covariance calculation */
    X_centered = X_opt.each_row() - arma::mean(X_opt, 0);
    Sigma_curr = (X_centered.t() * X_centered) / (n - 1.0);

    /* convergence check (Frobenius distance to target) */
    final_frob = arma::norm(Sigma_curr - Sigma_fixed, "fro");
    iter = i + 1;
    if (final_frob < tol) {
      converged = true;
      break;
    }

    /* covariance-constraint gradient (un-normalised) */
    grad_cov = 2.0 * X_centered * (Sigma_curr - Sigma_fixed);
    update_step = (lr * lambda) * grad_cov;

    /* gradient step (masked: only update originally-missing cells) */
    X_opt -= (update_step % mask);

    /* guard against numerical blow-up */
    if (X_opt.has_nan() || X_opt.has_inf()) {
      Rcpp::stop("Divergence detected: NaN or Inf produced during gradient descent.");
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("X_refined")  = X_opt,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("converged")  = converged,
    Rcpp::Named("final_frob") = final_frob
  );
}
