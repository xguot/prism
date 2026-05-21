library(smriti)
library(Matrix)

set.seed(42)

test_scaling <- function(n, p = 4, lambda = 0.1, lr = 0.001, max_iter = 1) {
  X <- matrix(rnorm(n * p), n, p)
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
  
  # Relative improvement
  rel_change <- (init_dist - final_dist) / init_dist
  
  return(list(n=n, init_dist=init_dist, final_dist=final_dist, change=change, rel_change=rel_change))
}

res_100 <- test_scaling(100)
res_10000 <- test_scaling(10000)

print(res_100)
print(res_10000)
