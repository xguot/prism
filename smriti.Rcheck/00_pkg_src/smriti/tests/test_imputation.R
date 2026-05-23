library(smriti)

# Test 1: Basic functionality and NA removal
cat("Running Test 1: Basic functionality and NA removal...\n")
df <- data.frame(A = c(1, 2, NA, 4), B = c(NA, 2, 3, 4))
result <- smriti::smriti_impute(df, time_cols = c("A", "B"), robust = FALSE)

# Assertion: All NAs must be gone
if (any(is.na(result))) {
  stop("Test 1 Failed: Imputed dataset still contains NAs.")
}
cat("Test 1 Passed.\n")

# Test 2: Structural Integrity (Observed values untouched)
cat("Running Test 2: Observed value preservation...\n")
# If observed values changed, crash the test
if (result$A[1] != 1) {
  stop(sprintf("Test 2 Failed: Observed value A[1] changed from 1 to %f", result$A[1]))
}
if (result$B[2] != 2) {
  stop(sprintf("Test 2 Failed: Observed value B[2] changed from 2 to %f", result$B[2]))
}
cat("Test 2 Passed.\n")

# Test 3: Custom Target implementation
cat("Running Test 3: Custom Target...\n")
target <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
result_custom <- smriti::smriti_impute(df, time_cols = c("A", "B"), custom_target = target)
if (any(is.na(result_custom))) {
  stop("Test 3 Failed: Custom target imputation contains NAs.")
}
cat("Test 3 Passed.\n")

cat("\nAll critical unit tests passed successfully.\n")
