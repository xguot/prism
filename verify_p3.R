library(smriti)

# Function to create a non-PSD Spearman correlation matrix from data
# We can use the fact that pairwise correlations can be inconsistent.
set.seed(42)
n <- 20
p <- 4

# Create data that is "almost" inconsistent
# V1 and V2 highly correlated
# V2 and V3 highly correlated
# V1 and V3 uncorrelated
# But we need missingness to make pairwise Spearman non-PSD

generate_non_psd_data <- function() {
  for (i in 1:1000) {
    X <- matrix(rnorm(n * p), n, p)
    # Introduce dependencies
    X[, 2] <- X[, 1] + rnorm(n, sd = 0.1)
    X[, 3] <- X[, 2] + rnorm(n, sd = 0.1)
    # X[, 1] and X[, 3] should be highly correlated, but we will force them not to be
    # by using missingness.
    
    # Randomly delete values
    X_miss <- X
    X_miss[sample(1:(n * p), 0.5 * n * p)] <- NA
    
    cor_mat <- cor(X_miss, method = "spearman", use = "pairwise.complete.obs")
    if (any(eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values < -1e-7)) {
      return(X_miss)
    }
  }
  stop("Failed to generate non-PSD data")
}

X_miss <- generate_non_psd_data()
df_miss <- as.data.frame(X_miss)
colnames(df_miss) <- paste0("V", 1:p)

cor_mat <- cor(df_miss, method = "spearman", use = "pairwise.complete.obs")
cat("Original Spearman correlation eigenvalues:\n")
print(eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values)

# Pass through robust=TRUE path
cat("\nRunning smriti_impute(robust = TRUE)...\n")
res <- try(smriti_impute(df_miss, time_cols = 1:p, robust = TRUE))

if (inherits(res, "try-error")) {
    cat("FAILED: smriti_impute crashed\n")
} else {
    cat("SUCCESS: smriti_impute completed\n")
    
    # Check if the target was actually projected
    # We can't easily check the internal target, but we can check if the result exists.
}
