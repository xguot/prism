library(prism)
library(lavaan)

# Test Data Generation
set.seed(42)
df_miss <- data.frame(
  T1 = c(1, 2, NA, 4, 5),
  T2 = c(NA, 2, 3, 4, 6),
  T3 = c(1, NA, 3, 4, 7)
)

# Define a simple growth model
model <- "
  i =~ 1*T1 + 1*T2 + 1*T3
  s =~ 0*T1 + 1*T2 + 2*T3
"

cat("Running Test 1: Basic FIML projection...\n")
result <- prism_fiml(df_miss, model)

if (any(is.na(result))) {
  stop("Test 1 Failed: Imputed dataset contains NAs.")
}
if (nrow(result) != nrow(df_miss) || ncol(result) != ncol(df_miss)) {
  stop("Test 1 Failed: Output dimensions do not match input.")
}
# Observed values must be untouched
obs_idx_T1 <- !is.na(df_miss$T1)
if (max(abs(result$T1[obs_idx_T1] - df_miss$T1[obs_idx_T1])) > 1e-12) {
  stop("Test 1 Failed: Observed values were modified.")
}
cat("Test 1 Passed.\n")

cat("Running Test 2: FIML projection with initial imputation...\n")
initial <- df_miss
for (j in 1:3) {
  initial[is.na(initial[, j]), j] <- mean(initial[, j], na.rm = TRUE)
}
result2 <- prism_fiml(df_miss, model, initial_imputation = initial)
if (any(is.na(result2))) {
  stop("Test 2 Failed: Imputed dataset contains NAs.")
}
cat("Test 2 Passed.\n")

cat("Running Test 3: Convergence with default lambda...\n")
result3 <- prism_fiml(df_miss, model, lambda = 1.0, tol = 1e-4)
if (any(is.na(result3))) {
  stop("Test 3 Failed: Imputed dataset contains NAs.")
}
cat("Test 3 Passed.\n")

cat("Running Test 4: Numerical guard — nearest PSD projection...\n")
non_psd <- matrix(c(1, 2, 2, 1), 2, 2)
psd_fixed <- prism:::nearest_psd(non_psd)
eig_vals  <- eigen(psd_fixed, only.values = TRUE)$values
if (any(eig_vals < -1e-12)) {
  stop("Test 4 Failed: nearest_psd did not yield a positive semidefinite matrix.")
}
cat("Test 4 Passed.\n")

cat("Running Test 5: Error on 100% missing column...\n")
df_broken <- df_miss
df_broken$T1 <- as.numeric(NA)
err_msg <- tryCatch(
  prism_fiml(df_broken, model),
  error = function(e) e$message
)
if (!grepl("100% missing", err_msg)) {
  stop("Test 5 Failed: Did not catch 100% missing column error. Got: ", err_msg)
}
cat("Test 5 Passed.\n")

cat("Running Test 6: Error on non-numeric column...\n")
df_bad <- df_miss
df_bad$T1 <- as.character(df_bad$T1)
err_msg2 <- tryCatch(
  prism_fiml(df_bad, model),
  error = function(e) e$message
)
if (!grepl("numeric", err_msg2)) {
  stop("Test 6 Failed: Did not catch non-numeric column. Got: ", err_msg2)
}
cat("Test 6 Passed.\n")

cat("Running Test 7: Multiple imputation via parameter perturbation...\n")
fit <- lavaan::growth(model, data = df_miss, missing = "fiml")
mi_list <- prism_mi(df_miss, fit, m = 3)

if (!inherits(mi_list, "prism_mi_list")) {
  stop("Test 7 Failed: Output is not of class 'prism_mi_list'.")
}
if (length(mi_list) != 3) {
  stop(sprintf("Test 7 Failed: Expected 3 imputations, got %d.", length(mi_list)))
}
if (any(sapply(mi_list, function(x) any(is.na(x))))) {
  stop("Test 7 Failed: One or more MI datasets contain NAs.")
}
# Check that the imputations are not identical (parameter perturbation worked)
cov1 <- stats::cov(mi_list[[1]][, c("T1", "T2", "T3")])
cov2 <- stats::cov(mi_list[[2]][, c("T1", "T2", "T3")])
if (isTRUE(all.equal(cov1, cov2, tolerance = 1e-8))) {
  stop("Test 7 Failed: Multiple imputations produced identical covariance matrices.")
}
cat("Test 7 Passed.\n")

cat("\nAll professional unit tests passed successfully.\n")
