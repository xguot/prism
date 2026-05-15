library(smriti)
library(missForest)
library(lavaan)

# Structural Validation Pipeline for ADNI Clinical Data
# Purpose: Evaluate variance recovery on empirical longitudinal trajectories.

# Path to the clinical dataset (assumed to be pre-processed and clean)
data_path <- "data/adni_clean.csv"

if (!file.exists(data_path)) {
        stop(paste("Clinical data not found at:", data_path, ". Ensure data is synced to Rivanna scratch."))
}

# Load data and isolate complete cases to establish the 'True' structural manifold
full_data <- read.csv(data_path)
complete_data <- full_data[complete.cases(full_data), ]

if (nrow(complete_data) < 100) {
        stop("Insufficient complete cases for structural validation.")
}

# Define the longitudinal time points (assuming T1, T2, T3, T4 columns)
time_cols <- grep("^T[1-4]", colnames(complete_data), value = TRUE)
if (length(time_cols) < 4) {
        stop("Required longitudinal columns (T1-T4) not found in dataset.")
}

# Calculate the 'True' population variance from the complete cases
# We use this as the benchmark for recovery after artificial corruption.
true_cov <- stats::cov(complete_data[, time_cols])
true_slope_var <- true_cov[2, 2] # Simplified proxy for demonstration

cat(sprintf("Baseline Complete Cases: %d\n", nrow(complete_data)))
cat(sprintf("Target Slope Variance (Complete): %.4f\n", true_slope_var))

# Apply a 30% MCAR mask to corrupt the complete cases
set.seed(42)
n_cells <- length(as.matrix(complete_data[, time_cols]))
mask_indices <- sample(1:n_cells, size = floor(0.3 * n_cells))
corrupted_matrix <- as.matrix(complete_data[, time_cols])
corrupted_matrix[mask_indices] <- NA
corrupted_data <- complete_data
corrupted_data[, time_cols] <- corrupted_matrix

# 1. missForest Imputation (ML baseline)
cat("Running missForest imputation...\n")
imp_mf <- missForest(corrupted_data)$ximp
mf_cov <- stats::cov(imp_mf[, time_cols])
mf_slope_var <- mf_cov[2, 2]

# 2. Smriti Imputation (Lagrangian Refinement)
# We use robust = FALSE here to prioritize efficiency on the established manifold
cat("Running smriti refinement (robust = FALSE)...\n")
imp_sm <- smriti_impute(corrupted_data, time_cols = time_cols, robust = FALSE)
sm_cov <- stats::cov(imp_sm[, time_cols])
sm_slope_var <- sm_cov[2, 2]

# Output Results
results <- data.frame(
        Method = c("Truth", "missForest", "Smriti"),
        Slope_Var = c(true_slope_var, mf_slope_var, sm_slope_var),
        Bias = c(0, mf_slope_var - true_slope_var, sm_slope_var - true_slope_var)
)

print(results)
write.csv(results, "empirical/adni_validation_results.csv", row.names = FALSE)
cat("Results saved to empirical/adni_validation_results.csv\n")
