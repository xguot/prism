library(smriti)
library(missForest)
library(lavaan)

# Structural Validation Pipeline for ADNI Clinical Data
# Refactored for Local Execution with Z-Score Standardization

# Path to the clinical dataset
data_path <- "data/adni_clean.csv"

if (!file.exists(data_path)) {
        stop(paste("Clinical data not found at:", data_path))
}

# Load data and isolate complete cases
full_data <- read.csv(data_path)
complete_data <- full_data[complete.cases(full_data), ]

# Define longitudinal columns
time_cols <- grep("^T[1-4]", colnames(complete_data), value = TRUE)

# Calculate True Variance from raw complete cases
true_cov <- stats::cov(complete_data[, time_cols])
true_slope_var <- true_cov[2, 2]

# Z-Score Normalization to prevent gradient explosion
# We scale the longitudinal columns to mean 0, variance 1
scaled_obj <- scale(complete_data[, time_cols])
center_attr <- attr(scaled_obj, "scaled:center")
scale_attr <- attr(scaled_obj, "scaled:scale")

complete_scaled <- complete_data
complete_scaled[, time_cols] <- scaled_obj

# Apply 30% MCAR mask to the scaled data
set.seed(42)
n_cells <- length(as.matrix(complete_scaled[, time_cols]))
mask_indices <- sample(1:n_cells, size = floor(0.3 * n_cells))
corrupted_matrix <- as.matrix(complete_scaled[, time_cols])
corrupted_matrix[mask_indices] <- NA
corrupted_scaled <- complete_scaled
corrupted_scaled[, time_cols] <- corrupted_matrix

# 1. missForest Imputation
cat("Running missForest imputation on scaled data...\n")
imp_mf_scaled <- missForest(corrupted_scaled)$ximp

# 2. Smriti Imputation (Refinement)
cat("Running smriti refinement on scaled data...\n")
# Resetting lambda to a more aggressive 0.5 since data is now unit-variance
imp_sm_scaled <- smriti_impute(corrupted_scaled, time_cols = time_cols, robust = FALSE, lambda = 0.5)

# Un-scale results back to raw magnitude (~7000)
unscale <- function(scaled_df, center, scale_val) {
        raw_mat <- t(apply(as.matrix(scaled_df[, time_cols]), 1, function(x) x * scale_val + center))
        colnames(raw_mat) <- time_cols
        res <- scaled_df
        res[, time_cols] <- raw_mat
        return(res)
}

imp_mf_raw <- unscale(imp_mf_scaled, center_attr, scale_attr)
imp_sm_raw <- unscale(imp_sm_scaled, center_attr, scale_attr)

# Calculate Final Variances
mf_slope_var <- stats::cov(imp_mf_raw[, time_cols])[2, 2]
sm_slope_var <- stats::cov(imp_sm_raw[, time_cols])[2, 2]

# Output Results
results <- data.frame(
        Method = c("Truth", "missForest", "Smriti"),
        Slope_Var = c(true_slope_var, mf_slope_var, sm_slope_var),
        Bias = c(0, mf_slope_var - true_slope_var, sm_slope_var - true_slope_var)
)

print(results)
write.csv(results, "empirical/adni_validation_results.csv", row.names = FALSE)
cat("Results saved to empirical/adni_validation_results.csv\n")
