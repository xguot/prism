#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Project the non-parametric imputation estimates back toward the structural
 * manifold defined by the target covariance. This Lagrangian constraint
 * prevents the variance collapse typically observed in random forest
 * imputations by penalizing deviations from the established longitudinal
 * variance structure.
 */
// [[Rcpp::export]]
arma::mat constrain_covariance(arma::mat X_imp, arma::mat Sigma_target, double lambda,
                         double lr, int max_iter)
{
  int n = X_imp.n_rows;
  arma::mat X_opt = X_imp;
  arma::mat Sigma_curr;
  arma::mat grad;

  for (int i = 0; i < max_iter; i++) {
    Sigma_curr = arma::cov(X_opt);

    /*
     * Normalize the gradient by the degrees of freedom (n-1) to ensure
     * the optimization step is invariant to sample size, stabilizing
     * the manifold projection across varying dataset scales.
     */
    grad = (X_opt * (Sigma_curr - Sigma_target)) / (n - 1);

    X_opt = X_opt - (lr * lambda) * grad;

    /*
     * Halt iteration once the Frobenius norm of the covariance deviation
     * falls below the precision threshold, indicating that the imputed
     * data has successfully converged onto the target manifold.
     */
    if (arma::norm(Sigma_curr - Sigma_target, "fro") < 1e-4)
      break;
  }

  return X_opt;
}
