library(smriti)

set.seed(42)
n <- 100
p <- 4
X <- matrix(rnorm(n * p), n, p)
X_miss <- X
X_miss[sample(1:(n * p), 0.2 * n * p)] <- NA

df_miss <- as.data.frame(X_miss)
colnames(df_miss) <- paste0("V", 1:p)

# Run smriti_impute
res <- smriti_impute(df_miss, time_cols = 1:p, robust = FALSE)

# Check observed-data invariance
obs_mask <- !is.na(X_miss)
drift <- abs(as.matrix(res[, 1:p])[obs_mask] - X[obs_mask])
max_drift <- max(drift)

cat(sprintf("Max drift in observed values: %e\n", max_drift))

if (max_drift > 0) {
    idx <- which(drift == max_drift, arr.ind = TRUE)
    # We need to find the coordinates in the original matrix
    obs_indices <- which(obs_mask, arr.ind = TRUE)
    target_idx <- obs_indices[idx[1], ]
    cat(sprintf("Drift at row %d, col %d: %e\n", target_idx[1], target_idx[2], max_drift))
}

# Bitwise identical check
is_identical <- identical(as.matrix(res[, 1:p])[obs_mask], X[obs_mask])
cat(sprintf("Bitwise identical: %s\n", is_identical))
