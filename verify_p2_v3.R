library(smriti)
library(Matrix)

set.seed(42)

test_scaling_fixed_dev <- function(n, p = 4, lambda = 0.1, lr = 0.001, max_iter = 1) {
  # Start with data that has a fixed covariance (not identity)
  Sigma_start <- diag(p)
  Sigma_start[1, 2] <- Sigma_start[2, 1] <- 0.5
  
  # Generate data with this covariance
  L <- chol(Sigma_start)
  X <- matrix(rnorm(n * p), n, p) %*% L
  
  mask <- matrix(1, n, p)
  Sigma_target <- diag(p)
  
  Sigma_init <- cov(X)
  init_dist <- sqrt(sum((Sigma_init - Sigma_target)^2))
  
  X_refined <- smriti:::constrain_covariance(
    X_imp = X,
    mask = mask,
    Sigma_target = Sigma_target,
    lambda = lambda,
    lr = lr,
    max_iter = max_iter
  )
  
  Sigma_final <- cov(X_refined)
  final_dist <- sqrt(sum((Sigma_final - Sigma_target)^2))
  
  change <- init_dist - final_dist
  
  return(list(n=n, init_dist=init_dist, final_dist=final_dist, change=change))
}

res_100 <- test_scaling_fixed_dev(100)
res_10000 <- test_scaling_fixed_dev(10000)

print(res_100)
print(res_10000)
