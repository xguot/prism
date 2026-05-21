library(smriti)
library(Matrix)

set.seed(42)

test_scaling <- function(n, p = 4, lambda = 0.1, lr = 0.001, max_iter = 1) {
  X <- matrix(rnorm(n * p), n, p)
  mask <- matrix(1, n, p) # All missing to allow movement
  Sigma_target <- diag(p)
  
  # Initial covariance
  Sigma_init <- cov(X)
  init_norm <- sqrt(sum((Sigma_init - Sigma_target)^2))
  
  # One iteration of gradient descent
  # We need to access the internal function. 
  # Since it's exported in NAMESPACE (indirectly via RcppExports), 
  # but maybe not from the package namespace if not in NAMESPACE file.
  # Let's check NAMESPACE.
  
  X_refined <- smriti:::constrain_covariance(
    X_imp = X,
    mask = mask,
    Sigma_target = Sigma_target,
    lambda = lambda,
    lr = lr,
    max_iter = max_iter
  )
  
  Sigma_final <- cov(X_refined)
  final_norm <- sqrt(sum((Sigma_final - Sigma_target)^2))
  
  change <- init_norm - final_norm
  return(change)
}

change_100 <- test_scaling(100)
change_10000 <- test_scaling(10000)

cat(sprintf("n=100, change in Frobenius norm: %e\n", change_100))
cat(sprintf("n=10000, change in Frobenius norm: %e\n", change_10000))

ratio <- change_10000 / change_100
cat(sprintf("Ratio (10000/100): %f\n", ratio))
