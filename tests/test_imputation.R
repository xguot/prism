library(smriti)

# Test Data Generation
set.seed(42)
df_miss <- data.frame(
  T1 = c(1, 2, NA, 4, 5),
  T2 = c(NA, 2, 3, 4, 6),
  T3 = c(1, NA, 3, 4, 7)
)

# Functional validation: Robust vs Non-Robust
cat("Running Test 1: Single Imputation (Robust vs Non-Robust)...\n")
res_nr <- smriti_impute(df_miss, time_cols = 1:3, robust = FALSE)
res_r  <- smriti_impute(df_miss, time_cols = 1:3, robust = TRUE)

if (any(is.na(res_nr)) || any(is.na(res_r))) {
  stop("Test 1 Failed: Imputed datasets contain NAs.")
}
if (identical(res_nr, res_r)) {
  stop("Test 1 Failed: Robust and Non-Robust paths yielded identical results on noisy data.")
}
cat("Test 1 Passed.\n")

# Multiple Imputation (MI) list structure and bootstrap integrity
cat("Running Test 2: Multiple Imputation (smriti_mi)...\n")
m <- 3
mi_list <- smriti_mi(df_miss, time_cols = 1:3, m = m, robust = TRUE)

if (!inherits(mi_list, "smriti_mi_list")) {
  stop("Test 2 Failed: Output is not of class 'smriti_mi_list'.")
}
if (length(mi_list) != m) {
  stop(sprintf("Test 2 Failed: Expected %d imputations, got %d.", m, length(mi_list)))
}
if (any(sapply(mi_list, function(x) any(is.na(x))))) {
  stop("Test 2 Failed: One or more MI datasets contain NAs.")
}
cat("Test 2 Passed.\n")

# Wrapper stability: missForest refinement logic
cat("Running Test 3: Wrapper Integrity (smriti_forest)...\n")
# Note: smriti_forest handles its own fallback if missForest is missing
res_forest <- smriti_forest(df_miss, time_cols = 1:3, robust = TRUE)
if (any(is.na(res_forest))) {
  stop("Test 3 Failed: smriti_forest result contains NAs.")
}
cat("Test 3 Passed.\n")

# Numerical guard: Nearest PSD projection for indefinite matrices
cat("Running Test 4: Nearest PSD Projection...\n")
# Create a non-PSD matrix (negative eigenvalue)
non_psd <- matrix(c(1, 2, 2, 1), 2, 2)
# eigen(non_psd)$values are 3 and -1
psd_fixed <- smriti:::nearest_psd(non_psd)
eig_vals  <- eigen(psd_fixed, only.values = TRUE)$values
if (any(eig_vals < -1e-12)) {
  stop("Test 4 Failed: nearest_psd did not yield a positive semidefinite matrix.")
}
cat("Test 4 Passed.\n")

# Fallback mechanism: Mean imputation when initial_imputation is NULL
cat("Running Test 5: Fallback Mean Imputation...\n")
# Expect a warning when initial_imputation is NULL
suppressWarnings({
  res_mean <- smriti_impute(df_miss, time_cols = 1:3, initial_imputation = NULL)
})
if (any(is.na(res_mean))) {
  stop("Test 5 Failed: Fallback mean imputation failed to remove NAs.")
}
cat("Test 5 Passed.\n")

# Structural guard: Detection of 100% missing longitudinal columns
cat("Running Test 6: Guard for 100% missing columns...\n")
df_broken <- df_miss
df_broken$T1 <- as.numeric(NA)
# Use column names to ensure we target the right columns
err_msg <- tryCatch(smriti_impute(df_broken, time_cols = c("T1", "T2", "T3")), error = function(e) e$message)

if (!grepl("100% missing", err_msg)) {
  stop("Test 6 Failed: Did not catch 100% missing column error. Got: ", err_msg)
}
cat("Test 6 Passed.\n")

cat("\nAll professional unit tests passed successfully.\n")
