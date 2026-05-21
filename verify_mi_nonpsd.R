library(smriti)

# Generate data that is likely to produce non-PSD bootstrap samples
set.seed(42)
n <- 20
p <- 4
X <- matrix(rnorm(n * p), n, p)
X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
X[, 3] <- X[, 2] + rnorm(n, sd = 0.01)
X_miss <- X
X_miss[sample(1:(n * p), 0.6 * n * p)] <- NA
df_miss <- as.data.frame(X_miss)
colnames(df_miss) <- paste0("V", 1:p)

# Try smriti_mi
cat("Running smriti_mi...\n")
res <- try(smriti_mi(df_miss, time_cols = 1:p, m = 5))

if (inherits(res, "try-error")) {
    cat("smriti_mi FAILED:\n")
    print(res)
} else {
    cat("smriti_mi SUCCESS\n")
}
