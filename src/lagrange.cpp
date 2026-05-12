#include <Rcpp.h>

using namespace Rcpp;

/*
 * Constrain longitudinal covariance shift to prevent tree hallucination.
 */
// [[Rcpp::export]]
NumericVector optimize_lagrange(NumericVector x, double lambda) {
        int i;
        int n = x.size();
        NumericVector result(n);

        for (i = 0; i < n; i++)
                if (x[i] > lambda)
                        result[i] = x[i] - lambda;
                else
                        result[i] = 0;

        return result;
}
