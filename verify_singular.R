library(smriti)

set.seed(42)
n <- 100
p <- 4
X <- matrix(rnorm(n * p), n, p)
df <- as.data.frame(X)
colnames(df) <- paste0("V", 1:p)

# Create a singular target covariance
# We can do this by setting robust=FALSE and providing data with identical columns
df_singular <- df
df_singular$V2 <- df_singular$V1
df_singular$V1[1:50] <- NA
df_singular$V2[51:100] <- NA
# Now V1 and V2 have no common observed values? 
# No, that's not what I want.

# Actually, smriti_impute calculates the target itself.
# Let's just pass a singular initial imputation and see if it works.
# Or better, use robust=TRUE on data that will yield a singular matrix.

df_zero <- df
df_zero$V1 <- 0
df_zero$V2 <- 0
# Spearman correlation will be NA if variance is zero.
# stats::mad will be 0.

res <- try(smriti_impute(df_zero, time_cols = 1:p, robust = TRUE))
if (inherits(res, "try-error")) {
    cat("FAILED on zero variance\n")
} else {
    cat("SUCCESS on zero variance\n")
}
