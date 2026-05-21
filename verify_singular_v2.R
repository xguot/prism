library(smriti)

set.seed(42)
n <- 100
p <- 4
X <- matrix(rnorm(n * p), n, p)
# Make V1 and V2 identical to ensure singularity
X[, 2] <- X[, 1]
df <- as.data.frame(X)
colnames(df) <- paste0("V", 1:p)

# Run smriti_impute with robust=TRUE
# This should result in a singular sigma_target (rank 3)
res <- try(smriti_impute(df, time_cols = 1:p, robust = TRUE))

if (inherits(res, "try-error")) {
    cat("FAILED on singular target\n")
    print(res)
} else {
    cat("SUCCESS on singular target\n")
    # Check if target is actually singular
    sigma_target <- stats::cov(df[, 1:p], use = "pairwise.complete.obs")
    cat("Eigenvalues of target:\n")
    print(eigen(sigma_target)$values)
}
