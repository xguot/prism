library(smriti)
library(missForest)

# EHR Robustness Validation Pipeline
# Purpose: Benchmark Smriti's MCD-robust refinement against skewed EHR data.

# Path to the EHR proxy dataset
data_path <- "data/mimic_proxy.csv"

if (!file.exists(data_path)) {
        stop(paste("EHR proxy data not found at:", data_path))
}

# Load the skewed dataset
full_data <- read.csv(data_path)
time_cols <- grep("^HeartRate", colnames(full_data), value = TRUE)

# Calculate Baseline 'True' Variance (Note: This is biased by the 999s, 
# which is the point of the test)
true_var <- stats::var(as.vector(as.matrix(full_data[, time_cols])), na.rm = TRUE)

# Z-Score Normalization to prevent C++ Lagrangian gradient explosion
scaled_obj <- scale(full_data[, time_cols])
center_attr <- attr(scaled_obj, "scaled:center")
scale_attr <- attr(scaled_obj, "scaled:scale")

data_scaled <- full_data
data_scaled[, time_cols] <- scaled_obj

# Apply 30% MCAR mask
set.seed(42)
corrupted_scaled <- data_scaled
n_cells <- length(as.matrix(corrupted_scaled[, time_cols]))
mask_indices <- sample(1:n_cells, size = floor(0.3 * n_cells))
corrupted_matrix <- as.matrix(corrupted_scaled[, time_cols])
corrupted_matrix[mask_indices] <- NA
corrupted_scaled[, time_cols] <- corrupted_matrix

# 1. missForest Imputation (Standard Baseline)
cat("Running missForest on skewed EHR data...\n")
imp_mf_scaled <- missForest(corrupted_scaled)$ximp

# 2. Smriti Refinement (Robust = TRUE)
# The MCD estimator should ignore the 999-valued outliers when establishing the manifold.
cat("Running smriti refinement (robust = TRUE) on skewed EHR data...\n")
imp_sm_scaled <- smriti_impute(corrupted_scaled, time_cols = time_cols, robust = TRUE, lambda = 0.5)

# Un-scale results back to raw HeartRate magnitude
unscale <- function(scaled_df, center, scale_val, cols) {
        res <- scaled_df
        for (i in seq_along(cols)) {
                res[[cols[i]]] <- scaled_df[[cols[i]]] * scale_val[i] + center[i]
        }
        return(res)
}

imp_mf_raw <- unscale(imp_mf_scaled, center_attr, scale_attr, time_cols)
imp_sm_raw <- unscale(imp_sm_scaled, center_attr, scale_attr, time_cols)

# Calculate Recovered Variances
# We look at the variance of the trajectories to see if Smriti ignored the artifacts.
mf_var <- stats::var(as.vector(as.matrix(imp_mf_raw[, time_cols])))
sm_var <- stats::var(as.vector(as.matrix(imp_sm_raw[, time_cols])))

# Results Comparison
results <- data.frame(
        Metric = "HeartRate Variance",
        Skewed_Baseline = true_var,
        missForest = mf_var,
        Smriti_Robust = sm_var
)

print(results)
write.csv(results, "empirical/mimic_ehr_validation_results.csv", row.names = FALSE)
cat("Results saved to empirical/mimic_ehr_validation_results.csv\n")
