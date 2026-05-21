library(smriti)

set.seed(20260521)
n <- 200
p <- 4

sigma_true <- matrix(0.7, p, p)
diag(sigma_true) <- 1.0
X_full <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma_true)
colnames(X_full) <- paste0("V", 1:p)

missing_rate <- 0.10
n_cells <- n * p
n_missing <- floor(n_cells * missing_rate)
missing_idx <- sample(n_cells, n_missing)
X_missing <- X_full
X_missing[missing_idx] <- NA
df <- as.data.frame(X_missing)
colnames(df) <- paste0("V", 1:p)

cat("=== Check 3 & 5: Convergence Diagnostics + No-Drift Guarantee ===\n")
cat("Data:", n, "x", p, "with", n_missing, "missing cells (",
    round(100*n_missing/n_cells, 1), "%)\n\n")

captured_warnings <- character(0)
result <- withCallingHandlers(
  smriti_impute(data = df, time_cols = paste0("V", 1:p),
                max_iter = 2000, tol = 1e-6, robust = FALSE),
  warning = function(w) {
    captured_warnings <<- c(captured_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

# --- Check 3: Convergence Diagnostics ---
cat("--- Check 3: Convergence Warning ---\n")
cat("Number of warnings captured:", length(captured_warnings), "\n")

# Display all captured warnings
if (length(captured_warnings) > 0) {
  cat("All warnings:\n")
  for (i in seq_along(captured_warnings)) {
    cat("  [", i, "] ", captured_warnings[i], "\n", sep = "")
  }
}

convergence_warning <- grep("did not reach tolerance", captured_warnings, value = TRUE)
if (length(convergence_warning) > 0) {
  cat("\nConvergence warning found:\n")
  cat("  ", convergence_warning[1], "\n")
  
  check3a <- grepl("did not reach tolerance", convergence_warning[1])
  check3b <- grepl("Final Dist", convergence_warning[1])
  check3c <- grepl("Tol", convergence_warning[1])
  
  cat("  [", if(check3a) "PASS" else "FAIL", "] Warning states 'did not reach tolerance'\n")
  cat("  [", if(check3b) "PASS" else "FAIL", "] Warning includes 'Final Dist'\n")
  cat("  [", if(check3c) "PASS" else "FAIL", "] Warning includes 'Tol'\n")
  
  check3_pass <- all(check3a, check3b, check3c)
} else {
  cat("No convergence warning found. Checking actual convergence...\n")
  
  final_cov <- cov(as.matrix(result[, paste0("V", 1:p)]))
  cov_obs <- cov(as.matrix(df[, paste0("V", 1:p)]), use = "pairwise.complete.obs")
  target_eig <- eigen(cov_obs, symmetric = TRUE)
  vals <- target_eig$values
  vals[vals < 0] <- 0
  target_cov <- target_eig$vectors %*% diag(vals) %*% t(target_eig$vectors)
  target_cov <- (target_cov + t(target_cov)) / 2
  
  final_dist <- sqrt(sum((final_cov - target_cov)^2))
  cat("  Final distance to target:", format(final_dist, digits = 4), "(tol = 1e-6)\n")
  
  if (final_dist <= 1e-6) {
    cat("  [PASS] Convergence achieved within 2000 iterations - no warning needed\n")
    check3_pass <- TRUE
  } else {
    cat("  [FAIL] Did not converge AND no warning was issued\n")
    check3_pass <- FALSE
  }
}

cat("\n  Check 3 OVERALL: ", if(check3_pass) "PASS" else "FAIL", "\n\n", sep = "")

# --- Check 5: No-Drift Guarantee ---
cat("--- Check 5: Structural Integrity (No-Drift) ---\n")

X_raw <- as.matrix(df[, paste0("V", 1:p)])
X_imputed <- as.matrix(result[, paste0("V", 1:p)])

obs_before <- X_raw[!is.na(X_raw)]
obs_after  <- X_imputed[!is.na(X_raw)]

max_drift <- max(abs(obs_before - obs_after))
cat("  Number of observed cells:", length(obs_before), "\n")
cat("  Number of missing cells imputed:", sum(is.na(X_raw)), "\n")
cat("  Maximum observed-value drift:", format(max_drift, digits = 15), "\n")

missing_after <- X_imputed[is.na(X_raw)]
cat("  Missing cells imputed (non-NA in output):", all(!is.na(missing_after)), "\n")

if (max_drift == 0) {
  cat("  [PASS] Maximum drift is exactly 0 -- observed values strictly untouched\n")
  check5_pass <- TRUE
} else if (max_drift < 1e-12) {
  cat("  [PASS] Maximum drift < 1e-12 -- within floating-point tolerance\n")
  check5_pass <- TRUE
} else {
  cat("  [FAIL] Maximum drift exceeds acceptable tolerance:", max_drift, "\n")
  drift_positions <- which(abs(obs_before - obs_after) > 1e-12)
  cat("  Positions with drift:", length(drift_positions), "\n")
  check5_pass <- FALSE
}

cat("\n  Check 5 OVERALL: ", if(check5_pass) "PASS" else "FAIL", "\n\n", sep = "")

overall_pass <- check3_pass && check5_pass
quit(status = if(overall_pass) 0 else 1)
