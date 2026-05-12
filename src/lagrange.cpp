#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Project imputed data matrix onto the target covariance manifold.
 * This prevents the non-parametric Random Forest from hallucinating 
 * splits that collapse the established longitudinal variance.
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
                 * Compute gradient of the Lagrangian penalty:
                 * L = ||Sigma_curr - Sigma_target||^2_F
                 * Scaled by sample size to stabilize the update step.
                 */
                grad = (X_opt * (Sigma_curr - Sigma_target)) / (n - 1);

                X_opt = X_opt - (lr * lambda) * grad;

                if (arma::norm(Sigma_curr - Sigma_target, "fro") < 1e-4)
                        break;
        }

        return X_opt;
}
