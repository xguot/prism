# ============================================================================
# QA Validation Script for smriti
# Tests: 1 (100% Missing Crash), 3 (Convergence Warning), 5 (Data Integrity)
# ============================================================================

library(smriti)

cat("========================================\n")
cat("  smriti QA Validation Suite\n")
cat("========================================\n\n")

# ---- Test 1: 100% Missing Crash Fix ----
cat("--- Test 1: 100% Missing Column ---\n")
n <- 100
p <- 4
set.seed(42)
full_data <- matrix(rnorm(n * p), nrow = n, ncol = p)
full_data[, 2] <- NA  # column 2 is 100% missing
df <- as.data.frame(full_data)
names(df) <- paste0("V", 1:p)

result1 <- tryCatch(
  smriti_impute(df, time_cols = paste0("V", 1:p)),
  error = function(e) {
    cat("Caught error: ", e$message, "\n")
    if (grepl("100% missing", e$message)) {
      cat("TEST 1: PASS — Correctly stopped with '100% missing' message.\n\n")
      return("PASS")
    } else {
      cat("TEST 1: FAIL — Error message does not mention '100% missing'.\n\n")
      return("FAIL")
    }
  }
)

if (!is.character(result1) || result1 != "PASS") {
  cat("TEST 1: FAIL — Did not stop as expected.\n\n")
}

# ---- Test 3: Convergence Diagnostic Logic ----
cat("--- Test 3: Convergence Warning ---\n")
n <- 200
p <- 4
set.seed(123)
X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)

# Introduce 10% MCAR missingness
set.seed(456)
na_mask <- matrix(runif(n * p) < 0.10, nrow = n, ncol = p)
X_na <- X
X_na[na_mask] <- NA
df_sim <- as.data.frame(X_na)
names(df_sim) <- paste0("V", 1:p)

warnings_captured <- character()
result3 <- withCallingHandlers(
  smriti_impute(df_sim, time_cols = paste0("V", 1:p),
                max_iter = 2000, tol = 1e-6, robust = FALSE),
  warning = function(w) {
    warnings_captured <<- c(warnings_captured, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("Captured warnings:", length(warnings_captured), "\n")
for (w in warnings_captured) {
  cat("  WARNING: ", w, "\n")
}

# Check for convergence warning
conv_warning <- grep("Covariance projection did not reach tolerance",
                     warnings_captured, value = TRUE)
if (length(conv_warning) > 0) {
  cat("\nConvergence warning found:\n  ", conv_warning, "\n")
  
  # Check that Final Dist and Tol are displayed
  if (grepl("Final Dist:", conv_warning) && grepl("Tol:", conv_warning)) {
    cat("TEST 3: PASS — Convergence warning issued with Final Dist and Tol.\n\n")
  } else {
    cat("TEST 3: PARTIAL — Warning found but missing Final Dist or Tol formatting.\n\n")
  }
} else {
  cat("\nTEST 3: NOTE — No convergence warning (final_dist <= tol).")
  cat(" This may be environment-dependent.\n\n")
}

# ---- Test 5: Data Integrity Regression Check ----
cat("--- Test 5: Data Integrity ---\n")
obs_before <- X_na[!is.na(X_na)]
obs_after  <- as.matrix(result3[, paste0("V", 1:p)])[!is.na(X_na)]
max_drift <- max(abs(obs_before - obs_after), na.rm = TRUE)
cat("Max observed-value drift:", max_drift, "\n")

if (max_drift == 0) {
  cat("TEST 5: PASS — Observed values unchanged.\n\n")
} else if (max_drift < 1e-12) {
  cat("TEST 5: PASS — Observed values unchanged (drift < 1e-12).\n\n")
} else {
  cat("TEST 5: FAIL — Observed values modified! Max drift:", max_drift, "\n\n")
}

# ---- Summary ----
cat("========================================\n")
cat("  Validation Complete\n")
cat("========================================\n")
