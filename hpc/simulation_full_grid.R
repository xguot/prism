library(MASS)
library(smriti)
library(parallel)

# Neutralize multi-threaded BLAS to prevent contention during parallel forking
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

set.seed(20250523)

grid_n    <- c(100, 200, 500, 1000, 5000)
grid_miss <- c(0.05, 0.15, 0.30)
grid_dist <- c("Normal", "Outlier", "Lognormal")
grid_mech <- c("MAR", "MNAR")
n_sims    <- 500
t_points  <- 4

mu_i <- 6.0
mu_s <- 2.0
v_i  <- 1.0
v_s  <- 1.0
c_is <- 0.0
v_e  <- 1.0

true_cov <- matrix(0, t_points, t_points)
for (r in 1:t_points) {
  for (c in 1:t_points) {
    true_cov[r, c] <- v_i + (r - 1) * (c - 1) * v_s + ((r - 1) + (c - 1)) * c_is
    if (r == c) true_cov[r, c] <- true_cov[r, c] + v_e
  }
}

frob_dist <- function(m1, m2) {
  sqrt(sum((m1 - m2)^2))
}

generate_data <- function(n, dist) {
  latent_vars <- mvrnorm(n, mu = c(mu_i, mu_s),
                         Sigma = matrix(c(v_i, c_is, c_is, v_s), 2, 2))

  if (dist == "Lognormal") {
    latent_vars <- exp(latent_vars)
    latent_vars[,1] <- scale(latent_vars[,1]) * sqrt(v_i) + mu_i
    latent_vars[,2] <- scale(latent_vars[,2]) * sqrt(v_s) + mu_s
  }

  data_mat <- matrix(0, n, t_points)
  for (j in 1:t_points) {
    err <- if(dist == "Lognormal") {
      scale(exp(rnorm(n))) * sqrt(v_e)
    } else {
      rnorm(n, 0, sqrt(v_e))
    }
    data_mat[, j] <- latent_vars[, 1] + (j - 1) * latent_vars[, 2] + err
  }

  if (dist == "Outlier") {
    out_idx <- sample(seq_len(n), floor(0.05 * n))
    data_mat[out_idx, ] <- data_mat[out_idx, ] + 5.0
  }

  df <- as.data.frame(data_mat)
  colnames(df) <- paste0("T", 1:t_points)
  df$true_slope <- latent_vars[, 2]

  return(df)
}

apply_missingness <- function(df, rate, mech) {
  df_miss <- df
  n <- nrow(df)

  if (mech == "MAR") {
    miss_n <- round(2 * n * rate / (t_points - 1))
    for (t in 1:(t_points - 1)) {
      order_idx <- order(df_miss[, t], decreasing = TRUE)
      target_rows <- (n - t * miss_n + 1):n
      if (target_rows[1] <= n) {
        real_idx <- order_idx[target_rows]
        df_miss[real_idx, (t + 1)] <- NA
      }
    }
  } else if (mech == "MNAR") {
    cor_ab <- 0.8
    a <- cor_ab / sqrt(1 - cor_ab^2)
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

run_iteration <- function(sim_id, params) {
  suppressPackageStartupMessages({
    library(MASS)
    library(lavaan)
    library(mice)
    library(ranger)
    library(smriti)
  })

  has_mf <- requireNamespace("missForest", quietly = TRUE)
  res_list <- list()

  df_true <- generate_data(params$n, params$dist)
  df_miss <- apply_missingness(df_true, params$miss, params$mech)

  file_prefix <- sprintf("sim_raw_data/cond_N%d_M%.2f_D%s_M%s_sim%03d", 
                         params$n, params$miss, params$dist, params$mech, sim_id)
  saveRDS(df_true, file = paste0(file_prefix, "_true.rds"))
  saveRDS(df_miss, file = paste0(file_prefix, "_miss.rds"))

  gcm_mod <- "
    i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
    s =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
    i ~~ s
    i ~~ i
    s ~~ s
  "

  time_fiml <- system.time({
    fit_fiml <- tryCatch(growth(gcm_mod, data = df_miss, missing = "fiml"), error = function(e) NULL)
    s_var_f <- NA; s_se_f <- NA; d_fiml <- NA
    if (!is.null(fit_fiml)) {
      pt <- parameterEstimates(fit_fiml)
      row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
      if(nrow(row) > 0) { s_var_f <- row$est[1]; s_se_f <- row$se[1] }
      imp_fiml <- tryCatch({
        pred_fiml <- lavaan::lavPredict(fit_fiml, type = "yhat")
        x_fiml <- as.matrix(df_miss[, 1:t_points])
        mask_fiml <- is.na(x_fiml)
        x_fiml[mask_fiml] <- pred_fiml[mask_fiml]
        x_fiml
      }, error = function(e) NULL)
      if (!is.null(imp_fiml)) d_fiml <- frob_dist(stats::cov(imp_fiml), true_cov)
    }
  })["elapsed"]

  res_list[[length(res_list) + 1]] <- data.frame(
    sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist, mech = params$mech,
    method = "FIML", f_dist = d_fiml, s_var = s_var_f, s_se = s_se_f, time_sec = unname(time_fiml)
  )

  time_mice <- system.time({
    imp_mice <- tryCatch({
      m_out <- mice::mice(df_miss, m = 1, method = "cart", printFlag = FALSE)
      mice::complete(m_out, 1)
    }, error = function(e) NULL)
    s_var_m <- NA; s_se_m <- NA; d_m <- NA
    if (!is.null(imp_mice)) {
      d_m <- frob_dist(stats::cov(imp_mice[,1:t_points]), true_cov)
      fit_m <- tryCatch(growth(gcm_mod, data = imp_mice), error = function(e) NULL)
      if (!is.null(fit_m)) {
        pt <- parameterEstimates(fit_m)
        row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_m <- row$est[1]; s_se_m <- row$se[1] }
      }
    }
  })["elapsed"]
  
  res_list[[length(res_list) + 1]] <- data.frame(
    sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist, mech = params$mech,
    method = "MICE", f_dist = d_m, s_var = s_var_m, s_se = s_se_m, time_sec = unname(time_mice)
  )

  time_mf <- system.time({
    imp_mf <- tryCatch({
      if (has_mf) {
        missForest::missForest(df_miss, verbose=FALSE)$ximp
      } else {
        x_init <- df_miss
        for(j in 1:4) x_init[is.na(x_init[,j]), j] <- mean(df_miss[,j], na.rm=TRUE)
        imp_r_base <- x_init
        for(j in 1:4) {
          na_idx <- is.na(df_miss[,j])
          if(any(na_idx)) {
            mod <- ranger::ranger(dependent.variable.name = colnames(df_miss)[j],
                                  data = x_init[!na_idx, ], num.trees = 50, verbose = FALSE)
            imp_r_base[na_idx, j] <- predict(mod, x_init[na_idx, ])$predictions
          }
        }
        imp_r_base
      }
    }, error = function(e) NULL)
    s_var_mf <- NA; s_se_mf <- NA; d_mf <- NA
    if (!is.null(imp_mf)) {
      d_mf <- frob_dist(stats::cov(imp_mf[,1:t_points]), true_cov)
      fit_mf <- tryCatch(growth(gcm_mod, data = imp_mf), error = function(e) NULL)
      if (!is.null(fit_mf)) {
        pt <- parameterEstimates(fit_mf)
        row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_mf <- row$est[1]; s_se_mf <- row$se[1] }
      }
    }
  })["elapsed"]

  res_list[[length(res_list) + 1]] <- data.frame(
    sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist, mech = params$mech,
    method = "missForest", f_dist = d_mf, s_var = s_var_mf, s_se = s_se_mf, time_sec = unname(time_mf)
  )

  time_snr <- system.time({
    imp_snr <- tryCatch(smriti_impute(df_miss, time_cols = 1:t_points, initial_imputation = imp_mf, robust = FALSE), error = function(e) NULL)
    s_var_snr <- NA; s_se_snr <- NA; d_snr <- NA
    if (!is.null(imp_snr)) {
      d_snr <- frob_dist(stats::cov(imp_snr[,1:t_points]), true_cov)
      fit_snr <- tryCatch(growth(gcm_mod, data = imp_snr), error = function(e) NULL)
      if (!is.null(fit_snr)) {
        pt <- parameterEstimates(fit_snr)
        row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_snr <- row$est[1]; s_se_snr <- row$se[1] }
      }
    }
  })["elapsed"]
  
  res_list[[length(res_list) + 1]] <- data.frame(
    sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist, mech = params$mech,
    method = "Smriti_NR", f_dist = d_snr, s_var = s_var_snr, s_se = s_se_snr, time_sec = unname(time_snr)
  )

  time_sr <- system.time({
    imp_sr <- tryCatch(smriti_impute(df_miss, time_cols = 1:t_points, initial_imputation = imp_mf, robust = TRUE), error = function(e) NULL)
    s_var_sr <- NA; s_se_sr <- NA; d_sr <- NA
    if (!is.null(imp_sr)) {
      d_sr <- frob_dist(stats::cov(imp_sr[,1:t_points]), true_cov)
      fit_sr <- tryCatch(growth(gcm_mod, data = imp_sr), error = function(e) NULL)
      if (!is.null(fit_sr)) {
        pt <- parameterEstimates(fit_sr)
        row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
        if(nrow(row) > 0) { s_var_sr <- row$est[1]; s_se_sr <- row$se[1] }
      }
    }
  })["elapsed"]
  
  res_list[[length(res_list) + 1]] <- data.frame(
    sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist, mech = params$mech,
    method = "Smriti_Robust", f_dist = d_sr, s_var = s_var_sr, s_se = s_se_sr, time_sec = unname(time_sr)
  )

  rm(df_true, df_miss, imp_mice, imp_mf, imp_snr, imp_sr)
  do.call(rbind, res_list)
}

conditions <- expand.grid(n = grid_n, miss = grid_miss, dist = grid_dist, mech = grid_mech, stringsAsFactors = FALSE)
total_conditions <- nrow(conditions)

array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.na(array_id)) {
  chunk_size <- ceiling(total_conditions / 4)
  start_idx <- (array_id - 1) * chunk_size + 1
  end_idx   <- min(array_id * chunk_size, total_conditions)
  current_conditions <- conditions[start_idx:end_idx, , drop = FALSE]
  output_file <- sprintf("simulation_results_hpc_part%d.rds", array_id)
} else {
  current_conditions <- conditions
  output_file <- "simulation_results_hpc.rds"
}

# Determine core count: check SLURM allocation first, then cap at 10 for local 16GB RAM
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(slurm_cores)) slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
num_cores <- if (!is.na(slurm_cores)) slurm_cores else 10

cat(sprintf("Execution Mode: %s\n", if(is.na(array_id)) "Local/Standard" else "SLURM Array"))
cat(sprintf("Parallel Backend: mclapply (forking) across %d cores\n", num_cores))
cat(sprintf("Conditions: %d | Replications: %d\n\n", nrow(current_conditions), n_sims))

dir.create("sim_raw_data", showWarnings = FALSE)
all_results <- list()

for (i in seq_len(nrow(current_conditions))) {
  params <- current_conditions[i, ]
  cat(sprintf("[%s] Condition %d/%d: N=%d, Miss=%.2f, Dist=%s, Mech=%s\n",
              Sys.time(), i, nrow(current_conditions), params$n, params$miss, params$dist, params$mech))

  res_cond <- mclapply(seq_len(n_sims), function(sim_id) {
    run_iteration(sim_id, params)
  }, mc.cores = num_cores, mc.preschedule = TRUE)

  res_cond <- res_cond[sapply(res_cond, is.data.frame)]
  if (length(res_cond) > 0) {
    res_df <- do.call(rbind, res_cond)
    all_results[[i]] <- res_df
    saveRDS(do.call(rbind, all_results), output_file)
  }
  
  rm(res_cond)
  gc()
}

cat("\nSimulation complete. Results saved to:", output_file, "\n")
