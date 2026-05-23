/* 
 * Simulation Study: Tang & Tong (2025) Replication
 * 
 * This script evaluates the recovery of latent growth curve parameters and 
 * covariance structures under Missing At Random (MAR) mechanisms and 
 * contaminated normal (outlier) distributions.
 */

library(MASS)
library(smriti)

/* Load optional dependencies with safety checks */
suppressPackageStartupMessages({
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan required")
  if (!requireNamespace("mice", quietly = TRUE)) stop("mice required")
  if (!requireNamespace("missForest", quietly = TRUE)) {
    /* Fallback to ranger-based logic if randomForest/missForest compilation failed */
    message("missForest not available; using internal ranger fallback for smriti_forest")
  }
  library(lavaan)
})

/* --- Simulation Configuration --- */
set.seed(2025)
n_obs    <- 400
n_sims   <- 50
t_points <- 4

/* Ground Truth: Linear Growth Curve Parameters */
mu_i <- 6.0
mu_s <- 2.0
v_i  <- 1.0
v_s  <- 1.0
c_is <- 0.0
v_e  <- 1.0

/* Compute True Population Covariance Matrix (Sigma) */
/* Sigma_tu = Var(I) + (t-1)(u-1)Var(S) + Cov(I,S)(t+u-2) + Delta*Var(E) */
true_cov <- matrix(0, t_points, t_points)
for (r in 1:t_points) {
  for (c in 1:t_points) {
    true_cov[r, c] <- v_i + (r-1)*(c-1)*v_s + ((r-1)+(c-1))*c_is
    if (r == c) true_cov[r, c] <- true_cov[r, c] + v_e
  }
}

/* Helper: Introduce MAR missingness based on Wave 1 scores */
apply_mar <- function(df, rate = 0.20) {
  /* Use wave 1 (T1) as the predictor for missingness in T2-T4 */
  /* Higher T1 values increase probability of missingness */
  logit_p <- -2.0 + 0.5 * scale(df$T1)
  prob_m  <- 1 / (1 + exp(-logit_p))
  
  /* Target approximately 'rate' missingness across waves 2-4 */
  for (j in 2:t_points) {
    m_idx <- runif(nrow(df)) < (prob_m * rate / mean(prob_m))
    df[m_idx, j] <- NA
  }
  df
}

/* Helper: Calculate Frobenius Norm Distance */
frob_dist <- function(m1, m2) {
  sqrt(sum((m1 - m2)^2))
}

/* Results Accumulators */
results <- data.frame(
  condition = character(),
  method    = character(),
  f_dist    = numeric(),
  slope_var = numeric(),
  stringsAsFactors = FALSE
)

/* --- Monte Carlo Loop --- */
pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for (s in 1:n_sims) {
  
  /* 1. Generate Base Data (Normal) */
  latent_vars <- mvrnorm(n_obs, mu = c(mu_i, mu_s), 
                         Sigma = matrix(c(v_i, c_is, c_is, v_s), 2, 2))
  data_norm <- matrix(0, n_obs, t_points)
  for (j in 1:t_points) {
    data_norm[, j] <- latent_vars[, 1] + (j-1)*latent_vars[, 2] + rnorm(n_obs, 0, sqrt(v_e))
  }
  colnames(data_norm) <- paste0("T", 1:t_points)
  data_norm <- as.data.frame(data_norm)

  /* 2. Generate Outlier Data (5% Contaminated Normal) */
  data_out <- data_norm
  out_idx  <- sample(1:n_obs, floor(0.05 * n_obs))
  /* Add high-leverage structural noise to contaminated rows */
  data_out[out_idx, ] <- data_out[out_idx, ] + matrix(rnorm(length(out_idx)*t_points, 0, 10), 
                                                      nrow = length(out_idx))

  /* 3. Apply MAR Missingness */
  miss_norm <- apply_mar(data_norm)
  miss_out  <- apply_mar(data_out)

  /* 4. Model Definition for lavaan */
  gcm_mod <- "
    i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
    s =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
    i ~~ s
    i ~~ i
    s ~~ s
  "

  /* --- Estimator Comparison --- */
  
  /* A. FIML (Gold Standard for MAR + Normal) */
  fit_fiml <- tryCatch(
    growth(gcm_mod, data = miss_norm, missing = "fiml"),
    error = function(e) NULL
  )
  if (!is.null(fit_fiml)) {
    s_var <- parameterEstimates(fit_fiml)$est[parameterEstimates(fit_fiml)$lhs == "s" & 
                                              parameterEstimates(fit_fiml)$op == "~~" & 
                                              parameterEstimates(fit_fiml)$rhs == "s"]
    results <- rbind(results, data.frame(condition="Normal", method="FIML", 
                                         f_dist=NA, slope_var=s_var))
  }

  /* B. missForest Baseline */
  imp_mf <- tryCatch({
    if (requireNamespace("missForest", quietly = TRUE)) {
      missForest::missForest(miss_norm)$ximp
    } else {
      /* Internal smriti fallback if missForest is missing */
      smriti_ranger(miss_norm, time_cols = 1:4)
    }
  }, error = function(e) NULL)
  
  if (!is.null(imp_mf)) {
    d <- frob_dist(stats::cov(imp_mf[,1:4]), true_cov)
    fit_mf <- tryCatch(growth(gcm_mod, data = imp_mf), error = function(e) NULL)
    s_var <- if(!is.null(fit_mf)) parameterEstimates(fit_mf)$est[7] else NA
    results <- rbind(results, data.frame(condition="Normal", method="missForest", 
                                         f_dist=d, slope_var=s_var))
  }

  /* C. MICE Baseline (Single Imputation for simplicity in simulation) */
  imp_mice <- tryCatch({
    m_out <- mice::mice(miss_norm, m = 1, method = "pmm", printFlag = FALSE)
    mice::complete(m_out, 1)
  }, error = function(e) NULL)
  
  if (!is.null(imp_mice)) {
    d <- frob_dist(stats::cov(imp_mice[,1:4]), true_cov)
    fit_mice <- tryCatch(growth(gcm_mod, data = imp_mice), error = function(e) NULL)
    s_var <- if(!is.null(fit_mice)) parameterEstimates(fit_mice)$est[7] else NA
    results <- rbind(results, data.frame(condition="Normal", method="MICE", 
                                         f_dist=d, slope_var=s_var))
  }

  /* D. Smriti Non-Robust (Normal Data) */
  imp_snr <- tryCatch(smriti_forest(miss_norm, time_cols = 1:4, robust = FALSE), 
                      error = function(e) NULL)
  if (!is.null(imp_snr)) {
    d <- frob_dist(stats::cov(imp_snr[,1:4]), true_cov)
    fit_snr <- tryCatch(growth(gcm_mod, data = imp_snr), error = function(e) NULL)
    s_var <- if(!is.null(fit_snr)) parameterEstimates(fit_snr)$est[7] else NA
    results <- rbind(results, data.frame(condition="Normal", method="Smriti_NR", 
                                         f_dist=d, slope_var=s_var))
  }

  /* E. Smriti Robust (Outlier Data) */
  imp_sr <- tryCatch(smriti_forest(miss_out, time_cols = 1:4, robust = TRUE), 
                     error = function(e) NULL)
  if (!is.null(imp_sr)) {
    d <- frob_dist(stats::cov(imp_sr[,1:4]), true_cov)
    fit_sr <- tryCatch(growth(gcm_mod, data = imp_sr), error = function(e) NULL)
    s_var <- if(!is.null(fit_sr)) parameterEstimates(fit_sr)$est[7] else NA
    results <- rbind(results, data.frame(condition="Outlier", method="Smriti_R", 
                                         f_dist=d, slope_var=s_var))
  }

  setTxtProgressBar(pb, s)
}
close(pb)

/* --- Final Synthesis --- */
cat("\n/* Summary Statistics: Frobenius Distance to True Covariance */\n")
print(aggregate(f_dist ~ condition + method, data = results, FUN = mean))

cat("\n/* Summary Statistics: Latent Slope Variance (Truth = 1.0) */\n")
print(aggregate(slope_var ~ condition + method, data = results, FUN = mean))
