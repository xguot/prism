library(smriti)

# Generate simple test data with known covariance structure
set.seed(42)
n <- 100
true_cov <- matrix(c(4, 2, 1, 0.5,
                     2, 3, 1.5, 0.8,
                     1, 1.5, 2, 0.6,
                     0.5, 0.8, 0.6, 1), nrow = 4)
X <- MASS::mvrnorm(n, mu = c(10, 12, 14, 16), Sigma = true_cov)
df <- as.data.frame(X)
colnames(df) <- paste0("T", 1:4)

# Introduce 20% MCAR
set.seed(123)
df_miss <- df
for (j in 1:4) df_miss[sample(n, 20), j] <- NA

# Mean imputation as initial (simulating what missForest would produce)
x_init <- as.matrix(df_miss)
for (j in 1:4) {
  na_idx <- is.na(x_init[, j])
  x_init[na_idx, j] <- mean(x_init[, j], na.rm = TRUE)
}

cat("=== Test 1: robust=FALSE, lambda=1.0, lr=0.001 ===\n")
res1 <- smriti_impute(df_miss, time_cols = 1:4, robust = FALSE,
                      initial_imputation = x_init,
                      lambda = 1.0, learning_rate = 0.001, max_iter = 2000)
cat("Initial cov:\n"); print(round(cov(x_init), 3))
cat("Target cov:\n"); print(round(cov(df_miss[,1:4], use="pairwise"), 3))
cat("Refined cov:\n"); print(round(cov(as.matrix(res1[,1:4])), 3))
cat("True cov:\n"); print(round(true_cov, 3))

cat("\n=== Test 2: robust=TRUE ===\n")
res2 <- smriti_impute(df_miss, time_cols = 1:4, robust = TRUE,
                      initial_imputation = x_init,
                      lambda = 1.0, learning_rate = 0.001, max_iter = 2000)
cat("Refined cov (robust):\n"); print(round(cov(as.matrix(res2[,1:4])), 3))

cat("\n=== Test 3: smriti_mi ===\n")
mi_res <- smriti_mi(df_miss, time_cols = 1:4, m = 3,
                    initial_imputation = x_init,
                    lambda = 1.0, learning_rate = 0.001, max_iter = 500)
cat("Number of imputations:", length(mi_res), "\n")
cat("Class:", class(mi_res), "\n")

cat("\n=== Test 4: edge case - no initial imputation ===\n")
res3 <- smriti_impute(df_miss, time_cols = 1:4, robust = FALSE,
                      lambda = 1.0, learning_rate = 0.001, max_iter = 2000)
cat("Refined cov (mean-fallback):\n"); print(round(cov(as.matrix(res3[,1:4])), 3))

cat("\n=== All smoke tests passed ===\n")
