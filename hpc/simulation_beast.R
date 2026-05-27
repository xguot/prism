library(MASS)
library(smriti)
library(parallel)
library(lavaan)
library(mice)
library(ranger)
library(missForest)
library(missRanger)

# Neutralize multi-threaded BLAS to prevent resource contention
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# Beast Mode: 14 cores provides optimized utilization on 16GB RAM system (~1.1GB per worker)
num_cores <- min(14, parallel::detectCores() - 1)
set.seed(20250523)

# Simulation Grid Configuration
grid_n    <- c(100, 200, 500, 1000, 5000)
grid_miss <- c(0.05, 0.15, 0.30)
grid_dist <- c("Normal", "Outlier", "Lognormal")
grid_mech <- c("MAR", "MNAR")
n_sims    <- 500
t_points  <- 4

# True Population Parameters
mu_i <- 6.0; mu_s <- 2.0; v_i <- 1.0; v_s <- 1.0; c_is <- 0.0; v_e <- 1.0

# Compute True Population Covariance Matrix
true_cov <- matrix(0, t_points, t_points)
for (r in 1:t_points) {
  for (c in 1:t_points) {
    true_cov[r, c] <- v_i + (r - 1) * (c - 1) * v_s + ((r - 1) + (c - 1)) * c_is
    if (r == c) true_cov[r, c] <- true_cov[r, c] + v_e
  }
}

# Distance Metric: Frobenius Norm
frob_dist <- function(m1, m2) sqrt(sum((m1 - m2)^2))

# Data Generation Engine
generate_data <- function(n, dist) {
  latent_vars <- mvrnorm(n, mu = c(mu_i, mu_s), Sigma = matrix(c(v_i, c_is, c_is, v_s), 2, 2))
  
  if (dist == "Lognormal") {
    # Induce lognormality while attempting to preserve moments
    latent_vars <- exp(latent_vars)
    latent_vars[,1] <- scale(latent_vars[,1]) * sqrt(v_i) + mu_i
    latent_vars[,2] <- scale(latent_vars[,2]) * sqrt(v_s) + mu_s
  }
  
  data_mat <- matrix(0, n, t_points)
  for (j in 1:t_points) {
    err <- if(dist == "Lognormal") scale(exp(rnorm(n))) * sqrt(v_e) else rnorm(n, 0, sqrt(v_e))
    data_mat[, j] <- latent_vars[, 1] + (j - 1) * latent_vars[, 2] + err
  }
  
  if (dist == "Outlier") {
    # Add heavy-tailed noise to a subset of observations
    idx <- sample(seq_len(n), floor(0.05 * n))
    data_mat[idx, ] <- data_mat[idx, ] + 5.0
  }
  
  df <- as.data.frame(data_mat)
  colnames(df) <- paste0("T", 1:t_points)
  df$true_slope <- latent_vars[, 2]
  return(df)
}

# Missingness Engine
apply_missingness <- function(df, rate, mech) {
  df_miss <- df; n <- nrow(df)
  
  if (mech == "MAR") {
    # Probabilistic Logistic MAR (Safe against propagating NAs)
    for (t in 1:(t_points - 1)) {
      # Isolate rows where the current time point is observed
      idx <- which(!is.na(df_miss[, t]))
      if (length(idx) > 0) {
        x_prev <- scale(df_miss[idx, t])
        # Calibrate threshold to match target rate
        p_miss <- 1 / (1 + exp(-(x_prev - qnorm(1 - rate))))
        # Sample dropout only for valid rows
        drop_idx <- idx[rbinom(length(idx), 1, p_miss) == 1]
        df_miss[drop_idx, (t + 1)] <- NA
      }
    }
  } else if (mech == "MNAR") {
    # Latent-dependent Dropout
    cor_ab <- 0.8; a <- cor_ab / sqrt(1 - cor_ab^2)
    aux_var <- a * df$true_slope + rnorm(n, 0, 1)
    miss_rate_t <- 2 * rate / (t_points - 1)
    for(j in 2:t_points) {
      crit <- qnorm((1 - (j - 1) * miss_rate_t), mean = a * mu_s, sd = sqrt(a^2 + 1))
      df_miss[which(aux_var > crit), j] <- NA
    }
  }
  
  df_miss$true_slope <- NULL
  return(df_miss)
}

# Worker Function
run_iteration <- function(sim_id, params) {
  df_true <- generate_data(params$n, params$dist)
  df_miss <- apply_missingness(df_true, params$miss, params$mech)
  
  gcm_mod <- "i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4\ns =~ 0*T1 + 1*T2 + 2*T3 + 3*T4\ni ~~ s\ni ~~ i\ns ~~ s"
  res_list <- list()
  
  # Initialize imputation objects to handle failures gracefully
  imp_mice <- NULL; imp_mf <- NULL; imp_mr <- NULL; imp_sr <- NULL
  
  # Analysis: FIML Baseline
  time_fiml <- system.time({
    fit_fiml <- tryCatch(growth(gcm_mod, data = df_miss, missing = "fiml"), error = function(e) NULL)
    s_var_f <- NA; s_se_f <- NA; d_fiml <- NA
    if (!is.null(fit_fiml)) {
      pt <- parameterEstimates(fit_fiml)
      row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
      if(nrow(row) > 0) { s_var_f <- row$est[1]; s_se_f <- row$se[1] }
      
      # Extract model-implied covariance directly to avoid conditional expectation bias (lavPredict)
      implied_cov <- tryCatch(lavaan::lavInspect(fit_fiml, "cov.ov"), error = function(e) NULL)
      if (!is.null(implied_cov)) {
        d_fiml <- frob_dist(implied_cov, true_cov)
      }
    }
  })["elapsed"]
  res_list[[1]] <- data.frame(sim_id=sim_id, N=params$n, miss=params$miss, dist=params$dist, mech=params$mech,
                              method="FIML", f_dist=d_fiml, s_var=s_var_f, s_se=s_se_f, time_sec=unname(time_fiml))
  
  # Analysis: MICE Baseline
  time_mice <- system.time({
    imp_mice <- tryCatch(mice::complete(mice::mice(df_miss, m = 1, method = "cart", printFlag = FALSE), 1), error = function(e) NULL)
    s_var_m <- NA; s_se_m <- NA; d_m <- NA
    if (!is.null(imp_mice)) {
      d_m <- frob_dist(stats::cov(imp_mice[,1:4]), true_cov)
      fit_m <- tryCatch(growth(gcm_mod, data = imp_mice), error = function(e) NULL)
      if (!is.null(fit_m)) {
        pt <- parameterEstimates(fit_m); row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_m <- row$est[1]; s_se_m <- row$se[1] }
      }
    }
  })["elapsed"]
  res_list[[2]] <- data.frame(sim_id=sim_id, N=params$n, miss=params$miss, dist=params$dist, mech=params$mech,
                              method="MICE", f_dist=d_m, s_var=s_var_m, s_se=s_se_m, time_sec=unname(time_mice))

  # Analysis: missForest Baseline
  time_mf <- system.time({
    imp_mf <- tryCatch(missForest::missForest(df_miss, verbose=FALSE)$ximp, error = function(e) NULL)
    s_var_mf <- NA; s_se_mf <- NA; d_mf <- NA
    if (!is.null(imp_mf)) {
      d_mf <- frob_dist(stats::cov(imp_mf[,1:4]), true_cov)
      fit_mf <- tryCatch(growth(gcm_mod, data = imp_mf), error = function(e) NULL)
      if (!is.null(fit_mf)) {
        pt <- parameterEstimates(fit_mf); row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_mf <- row$est[1]; s_se_mf <- row$se[1] }
      }
    }
  })["elapsed"]
  res_list[[3]] <- data.frame(sim_id=sim_id, N=params$n, miss=params$miss, dist=params$dist, mech=params$mech,
                              method="missForest", f_dist=d_mf, s_var=s_var_mf, s_se=s_se_mf, time_sec=unname(time_mf))

  # Analysis: missRanger Baseline
  time_mr <- system.time({
    imp_mr <- tryCatch(missRanger::missRanger(df_miss, verbose=0), error = function(e) NULL)
    s_var_mr <- NA; s_se_mr <- NA; d_mr <- NA
    if (!is.null(imp_mr)) {
      d_mr <- frob_dist(stats::cov(imp_mr[,1:4]), true_cov)
      fit_mr <- tryCatch(growth(gcm_mod, data = imp_mr), error = function(e) NULL)
      if (!is.null(fit_mr)) {
        pt <- parameterEstimates(fit_mr); row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_mr <- row$est[1]; s_se_mr <- row$se[1] }
      }
    }
  })["elapsed"]
  res_list[[4]] <- data.frame(sim_id=sim_id, N=params$n, miss=params$miss, dist=params$dist, mech=params$mech,
                              method="missRanger", f_dist=d_mr, s_var=s_var_mr, s_se=s_se_mr, time_sec=unname(time_mr))
  
  # Analysis: Smriti Robust Projection
  time_sr <- system.time({
    imp_sr <- tryCatch(smriti_impute(df_miss, time_cols = 1:4, initial_imputation = imp_mf, robust = TRUE), error = function(e) NULL)
    s_var_sr <- NA; s_se_sr <- NA; d_sr <- NA
    if (!is.null(imp_sr)) {
      d_sr <- frob_dist(stats::cov(imp_sr[,1:4]), true_cov)
      fit_sr <- tryCatch(growth(gcm_mod, data = imp_sr), error = function(e) NULL)
      if (!is.null(fit_sr)) {
        pt <- parameterEstimates(fit_sr); row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_sr <- row$est[1]; s_se_sr <- row$se[1] }
      }
    }
  })["elapsed"]
  res_list[[5]] <- data.frame(sim_id=sim_id, N=params$n, miss=params$miss, dist=params$dist, mech=params$mech,
                              method="Smriti_Robust", f_dist=d_sr, s_var=s_var_sr, s_se=s_se_sr, time_sec=unname(time_sr))
  
  return(list(
    res = do.call(rbind, res_list), 
    raw_true = df_true, 
    raw_miss = df_miss,
    imp_mice = imp_mice,
    imp_mf   = imp_mf,
    imp_mr   = imp_mr,
    refined  = imp_sr
  ))
}

# Execution Pipeline
conditions <- expand.grid(n = grid_n, miss = grid_miss, dist = grid_dist, mech = grid_mech, stringsAsFactors = FALSE)

cat(sprintf("[%s] Standard Local Mode: Processing all %d conditions\n", as.character(Sys.time()), nrow(conditions)))

dir.create("sim_results", showWarnings = FALSE)
dir.create("sim_raw_data", showWarnings = FALSE)
dir.create("sim_plots", showWarnings = FALSE)

cat(sprintf("[%s] Launching Refined Local Mode across %d cores\n", as.character(Sys.time()), num_cores))

for (i in seq_len(nrow(conditions))) {
  params <- conditions[i, ]
  
  cat(sprintf("[%s] Condition %d/%d: N=%d, Miss=%.2f, Dist=%s, Mech=%s\n", 
              as.character(Sys.time()), i, nrow(conditions), params$n, params$miss, params$dist, params$mech))
  
  # Static chunking (mc.preschedule=TRUE) ensures stable Copy-On-Write for 16GB RAM
  out_list <- mclapply(seq_len(n_sims), function(s) run_iteration(s, params), mc.cores = num_cores, mc.preschedule = TRUE)
  
  # Data consolidation
  valid_out <- out_list[sapply(out_list, is.list)]
  if (length(valid_out) > 0) {
    results_df <- do.call(rbind, lapply(valid_out, `[[`, "res"))
    raw_data   <- lapply(valid_out, function(x) list(
      true     = x$raw_true, 
      miss     = x$raw_miss, 
      imp_mice = x$imp_mice, 
      imp_mf   = x$imp_mf, 
      imp_mr   = x$imp_mr,
      refined  = x$refined
    ))
    
    saveRDS(results_df, sprintf("sim_results/results_cond_%d.rds", i))
    saveRDS(raw_data,   sprintf("sim_raw_data/raw_data_cond_%d.rds", i))
    
    # Export Scatterplot Matrices for the first simulation of this condition
    pdf(sprintf("sim_plots/plots_cond_%d.pdf", i), width = 12, height = 12)
    par(mfrow = c(3, 2))
    for (name in c("true", "miss", "imp_mice", "imp_mf", "refined")) {
      d <- raw_data[[1]][[name]]
      if (!is.null(d)) {
        pairs(d[, 1:4], main = paste0("Condition ", i, ": ", name), col = rgb(0,0,0,0.2))
      }
    }
    dev.off()
    
    # Explicit Cleanup to stay under 16GB
    rm(results_df, raw_data, valid_out)
  }
  
  rm(out_list)
  gc()
}

cat("\nSimulation complete. All manuscript artifacts are analysis-ready.\n")
