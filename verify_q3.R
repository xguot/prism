library(smriti)
library(MASS)

set.seed(42)
n <- 200
p <- 4
# True covariance
Sigma_true <- matrix(0.5, p, p)
diag(Sigma_true) <- 1
# Generate data
data_clean <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma_true)
df <- as.data.frame(data_clean)
colnames(df) <- paste0("T", 1:p)

# Introduce 10% MCAR
mcar_idx <- sample(seq_len(n * p), 0.1 * n * p)
df_miss_mat <- as.matrix(df)
df_miss_mat[mcar_idx] <- NA
df_miss <- as.data.frame(df_miss_mat)

# Target covariance (from observed data)
target_cov <- cov(df_miss, use = "pairwise.complete.obs")
# Need to make sure it's PSD if we were doing it manually, but smriti_impute does it.

# Run smriti_impute
# lambda=0.1, learning_rate=0.001, max_iter=2000 are defaults
# but I'll specify them to be sure.
result <- smriti_impute(
  data = df_miss,
  time_cols = 1:p,
  lambda = 0.1,
  learning_rate = 0.001,
  max_iter = 2000,
  tol = 1e-6,
  robust = FALSE # Using Pearson for direct comparison with cov()
)

# Initial distance (smriti fallback is column mean)
x_raw <- as.matrix(df_miss)
x_init <- x_raw
for(i in 1:p) {
  x_init[is.na(x_init[,i]), i] <- mean(x_init[,i], na.rm = TRUE)
}
init_cov <- cov(x_init)
# smriti uses nearest_psd on target, so we should too for distance comparison
nearest_psd <- function(mat) {
  eig <- eigen(mat, symmetric = TRUE)
  vals <- eig$values
  vals[vals < 0] <- 0
  result <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  (result + t(result)) / 2
}
sigma_target <- nearest_psd(target_cov)

initial_dist <- sqrt(sum((init_cov - sigma_target)^2))
final_cov <- cov(as.matrix(result[, 1:p]))
final_dist <- sqrt(sum((final_cov - sigma_target)^2))

cat(sprintf("Initial Frobenius distance: %.6f\n", initial_dist))
cat(sprintf("Final Frobenius distance: %.6f\n", final_dist))

# Check observed values
obs_mask <- !is.na(as.matrix(df_miss))
obs_before <- as.matrix(df_miss)[obs_mask]
obs_after <- as.matrix(result[, 1:p])[obs_mask]
max_diff <- max(abs(obs_before - obs_after))
cat(sprintf("Max observed-value change: %.6e\n", max_diff))
