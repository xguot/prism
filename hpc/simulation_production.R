# Rivanna local library path (R_LIBS_USER may be ignored by site .Renviron)
.libPaths(c("~/R/rivanna-lib", .libPaths()))

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

# Core allocation for production run
num_cores <- min(24, parallel::detectCores() - 1)
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
seed_base <- 20250523
set.seed(if (is.na(array_id)) seed_base else seed_base + array_id)

# ── Simulation Grid ───────────────────────────────────────────────────────────
# Aligned with Tang & Tong (UVA) manuscript:
#   N       — 100, 200, 500, 1000, 5000, 10000
#   miss    — 5%, 10%, 15%, 30%
#   dist    — Normal, t(5), Outlier, Lognormal
#   mech    — MAR, MNAR
grid_n    <- c(100, 200, 500, 1000, 5000, 10000)
grid_miss <- c(0.05, 0.10, 0.15, 0.30)
grid_dist <- c("Normal", "t5", "Outlier", "Lognormal")
grid_mech <- c("MAR", "MNAR")
n_sims    <- 500
t_points  <- 4

# ── True Population Parameters ───────────────────────────────────────────────
#   σ²_i = 1, σ²_s = 1, σ_is = 0, σ²_e = 1, μ_i = 6, μ_s = 2
mu_i <- 6.0; mu_s <- 2.0; v_i <- 1.0; v_s <- 1.0; c_is <- 0.0; v_e <- 1.0

# True population covariance matrix (4 × 4)
true_cov <- matrix(0, t_points, t_points)
for (r in 1:t_points) {
  for (c in 1:t_points) {
    true_cov[r, c] <- v_i + (r - 1) * (c - 1) * v_s + ((r - 1) + (c - 1)) * c_is
    if (r == c) true_cov[r, c] <- true_cov[r, c] + v_e
  }
}

# ── Metrics ──────────────────────────────────────────────────────────────────
frob_dist <- function(m1, m2) sqrt(sum((m1 - m2)^2))
rel_bias  <- function(est, truth) {
  if (abs(truth) < 1e-12) return(est - truth)  # Tang & Tong: RB = raw bias when θ = 0
  100 * (est - truth) / truth
}

# ── Data Generation Engine ───────────────────────────────────────────────────
generate_data <- function(n, dist) {
  latent_vars <- mvrnorm(n, mu = c(mu_i, mu_s),
                         Sigma = matrix(c(v_i, c_is, c_is, v_s), 2, 2))

  data_mat <- matrix(0, n, t_points)
  for (j in 1:t_points) {
    err <- if (dist == "Lognormal") {
      scale(exp(rnorm(n))) * sqrt(v_e)
    } else if (dist == "t5") {
      # Student's t with 5 df, scaled to variance σ²_e = 1
      rt(n, df = 5) * sqrt(v_e * (5 - 2) / 5)
    } else {
      rnorm(n, 0, sqrt(v_e))
    }

    data_mat[, j] <- latent_vars[, 1] + (j - 1) * latent_vars[, 2] + err
  }

  if (dist == "Outlier") {
    # Subject-level measurement error outliers (5%)
    idx <- sample(seq_len(n), floor(0.05 * n))
    data_mat[idx, ] <- data_mat[idx, ] + 5.0
  }

  df <- as.data.frame(data_mat)
  colnames(df) <- paste0("T", 1:t_points)
  df$true_slope <- latent_vars[, 2]
  return(df)
}

# ── Missingness Engine ───────────────────────────────────────────────────────
apply_missingness <- function(df, rate, mech) {
  df_miss <- df; n <- nrow(df)

  if (mech == "MAR") {
    # Deterministic threshold MAR (Tang & Tong 2025 convention)
    # Low values on T_t cause missingness on T_{t+1}
    miss_n_per_step <- round(2 * n * rate / (t_points - 1))

    for (t in 1:(t_points - 1)) {
      # Identify subjects not already dropped at previous timepoints
      obs_idx <- which(!is.na(df_miss[, t]))
      if (length(obs_idx) > 0) {
        # Sort observed values: smallest values are at the end if decreasing=TRUE?
        # Actually, image says: order(..., decreasing=TRUE), target_rows <- (n - t*miss_n + 1):n
        # That targets the smallest values.
        order_idx <- order(df_miss[, t], decreasing = TRUE)
        # Calculate how many to drop this step (linear accumulation)
        target_rows <- (n - t * miss_n_per_step + 1):n
        target_rows <- target_rows[target_rows > 0 & target_rows <= n]

        if (length(target_rows) > 0) {
          drop_idx <- order_idx[target_rows]
          df_miss[drop_idx, (t + 1):t_points] <- NA # Dropout: once missing, stay missing
        }
      }
    }
  } else if (mech == "MNAR") {
    # Latent-slope-dependent dropout (Tang & Tong convention)
    cor_ab <- 0.8; a <- cor_ab / sqrt(1 - cor_ab^2)
    aux_var <- a * df$true_slope + rnorm(n, 0, 1)
    miss_rate_t <- 2 * rate / (t_points - 1)
    for (j in 2:t_points) {
      crit <- qnorm((1 - (j - 1) * miss_rate_t), mean = a * mu_s, sd = sqrt(a^2 + 1))
      df_miss[which(aux_var > crit), j] <- NA
    }
  }

  actual_miss <- mean(is.na(df_miss[, 1:t_points]))
  df_miss$true_slope <- NULL
  return(list(data = df_miss, actual_miss = actual_miss))
}

# ── Helper: pool multiple imputation results via Rubin's rules ────────────────
pool_rubin <- function(ests, ses) {
  m <- length(ests)
  if (m == 0) return(c(est = NA_real_, se = NA_real_))
  if (m == 1) return(c(est = ests[1], se = ses[1]))

  mean_est <- mean(ests, na.rm = TRUE)
  vw <- mean(ses^2, na.rm = TRUE)
  vb <- var(ests, na.rm = TRUE)
  total_var <- vw + (1 + 1/m) * vb
  c(est = mean_est, se = sqrt(total_var))
}

# ── Helper: extract GCM parameters from lavaan fit ────────────────────────────
# Returns a named vector of 5 structural parameters:
#   beta_L, beta_S, psi_L, psi_S, psi_LS
extract_gcm_params <- function(fit, type = "est") {
  pt <- parameterEstimates(fit)
  get_val <- function(lhs, op, rhs) {
    row <- pt[pt$lhs == lhs & pt$op == op & pt$rhs == rhs, ]
    if (nrow(row) > 0) row[[type]][1] else NA_real_
  }
  c(
    beta_L = get_val("i", "~1", ""),
    beta_S = get_val("s", "~1", ""),
    psi_L  = get_val("i", "~~", "i"),
    psi_S  = get_val("s", "~~", "s"),
    psi_LS = get_val("i", "~~", "s")
  )
}

# ── Helper: extract slope variance + SE ──────────────────────────────────────
extract_slope_var <- function(data, gcm_mod) {
  fit <- tryCatch(growth(gcm_mod, data = data), error = function(e) NULL)
  if (is.null(fit)) return(c(s_var = NA, s_se = NA))
  pt <- parameterEstimates(fit)
  row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
  if (nrow(row) > 0) c(s_var = row$est[1], s_se = row$se[1])
  else c(s_var = NA, s_se = NA)
}

# ── Single-Replication Worker ────────────────────────────────────────────────
run_iteration <- function(sim_id, params) {
  gcm_mod <- "
    i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
    s =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
    i ~~ s
    i ~~ i
    s ~~ s
  "

  make_row <- function(method, f_dist, s_var, s_se, time_sec,
                       beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA,
                       pipeline_time = time_sec) {
    data.frame(
      sim_id = sim_id, N = params$n, miss = params$miss,
      actual_miss = act_m,
      dist = params$dist, mech = params$mech,
      method    = method,
      f_dist    = f_dist,
      s_var     = s_var,
      s_var_bias = rel_bias(s_var, v_s),
      s_se      = s_se,
      est_L     = beta_L,   est_S     = beta_S,
      est_var_L = psi_L,    est_var_S = psi_S,    est_cov_LS = psi_LS,
      bias_L    = rel_bias(beta_L, mu_i),
      bias_S    = rel_bias(beta_S, mu_s),
      bias_var_L  = rel_bias(psi_L,  v_i),
      bias_var_S  = rel_bias(psi_S,  v_s),
      bias_cov_LS = rel_bias(psi_LS, c_is),
      time_sec     = time_sec,
      pipeline_time = pipeline_time,
      stringsAsFactors = FALSE
    )
  }
  res_list <- list()

  df_true <- generate_data(params$n, params$dist)
  miss_out <- apply_missingness(df_true, params$miss, params$mech)
  df_miss <- miss_out$data
  act_m   <- miss_out$actual_miss

  # ── FIML Baseline ────────────────────────────────────────────────────────
  time_fiml <- system.time({
    fit_fiml <- tryCatch(growth(gcm_mod, data = df_miss, missing = "fiml"),
                         error = function(e) NULL)
    s_var_f <- NA; s_se_f <- NA; d_fiml <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    if (!is.null(fit_fiml)) {
      pt <- parameterEstimates(fit_fiml)
      row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
      if (nrow(row) > 0) { s_var_f <- row$est[1]; s_se_f <- row$se[1] }
      implied_cov <- tryCatch(lavaan::lavInspect(fit_fiml, "cov.ov"),
                              error = function(e) NULL)
      if (!is.null(implied_cov)) d_fiml <- frob_dist(implied_cov, true_cov)
      gp <- extract_gcm_params(fit_fiml)
    }
  })["elapsed"]
  res_list[[1]] <- make_row("FIML", d_fiml, s_var_f, s_se_f, unname(time_fiml),
                            beta_L = gp["beta_L"], beta_S = gp["beta_S"],
                            psi_L = gp["psi_L"], psi_S = gp["psi_S"],
                            psi_LS = gp["psi_LS"])

  # ── MICE Baseline (MI m=20) ─────────────────────────────────────────────────
  time_mice <- system.time({
    m_mice <- 20
    imp_mice_list <- tryCatch({
      imp_obj <- mice::mice(df_miss, m = m_mice, method = "cart", printFlag = FALSE)
      mice::complete(imp_obj, "all")
    }, error = function(e) NULL)

    s_var_m <- NA; s_se_m <- NA; d_m <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)

    if (!is.null(imp_mice_list)) {
      # 1. Pool Covariance for Frobenius Distance
      cov_list <- lapply(imp_mice_list, function(x) stats::cov(x[, 1:t_points]))
      avg_cov <- Reduce("+", cov_list) / length(cov_list)
      d_m <- frob_dist(avg_cov, true_cov)

      # 2. Pool GCM Parameters using Rubin's Rules
      fit_results <- lapply(imp_mice_list, function(ds) {
        fit <- tryCatch(lavaan::growth(gcm_mod, data = ds), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        list(est = extract_gcm_params(fit, "est"), se = extract_gcm_params(fit, "se"))
      })

      valid_fits <- fit_results[!sapply(fit_results, is.null)]
      if (length(valid_fits) > 0) {
        p_names <- names(valid_fits[[1]]$est)
        pooled <- sapply(p_names, function(pn) {
          ests <- sapply(valid_fits, function(f) f$est[pn])
          ses  <- sapply(valid_fits, function(f) f$se[pn])
          pool_rubin(ests, ses)
        })
        gp <- pooled["est", ]
        gps <- pooled["se", ]
        s_var_m <- gp["psi_S"]; s_se_m <- gps["psi_S"]
      }
    }
  })["elapsed"]
  res_list[[2]] <- make_row("MICE", d_m, s_var_m, s_se_m, unname(time_mice),
                            beta_L = gp["beta_L"], beta_S = gp["beta_S"],
                            psi_L = gp["psi_L"], psi_S = gp["psi_S"],
                            psi_LS = gp["psi_LS"])

  # ── missForest Baseline ───────────────────────────────────────────────────
  time_mf <- system.time({
    imp_mf <- tryCatch(missForest::missForest(df_miss, verbose = FALSE)$ximp,
                       error = function(e) NULL)
    s_var_mf <- NA; s_se_mf <- NA; d_mf <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    if (!is.null(imp_mf)) {
      d_mf <- frob_dist(stats::cov(imp_mf[, 1:t_points]), true_cov)
      sv <- extract_slope_var(imp_mf, gcm_mod)
      s_var_mf <- sv["s_var"]; s_se_mf <- sv["s_se"]
      fit_mf <- tryCatch(growth(gcm_mod, data = imp_mf), error = function(e) NULL)
      if (!is.null(fit_mf)) gp <- extract_gcm_params(fit_mf)
    }
  })["elapsed"]
  res_list[[3]] <- make_row("missForest", d_mf, s_var_mf, s_se_mf, unname(time_mf),
                            beta_L = gp["beta_L"], beta_S = gp["beta_S"],
                            psi_L = gp["psi_L"], psi_S = gp["psi_S"],
                            psi_LS = gp["psi_LS"])

  # ── missRanger Baseline ───────────────────────────────────────────────────
  time_mr <- system.time({
    imp_mr <- tryCatch(missRanger::missRanger(df_miss, verbose = 0),
                       error = function(e) NULL)
    s_var_mr <- NA; s_se_mr <- NA; d_mr <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    if (!is.null(imp_mr)) {
      d_mr <- frob_dist(stats::cov(imp_mr[, 1:t_points]), true_cov)
      sv <- extract_slope_var(imp_mr, gcm_mod)
      s_var_mr <- sv["s_var"]; s_se_mr <- sv["s_se"]
      fit_mr <- tryCatch(growth(gcm_mod, data = imp_mr), error = function(e) NULL)
      if (!is.null(fit_mr)) gp <- extract_gcm_params(fit_mr)
    }
  })["elapsed"]
  res_list[[4]] <- make_row("missRanger", d_mr, s_var_mr, s_se_mr, unname(time_mr),
                            beta_L = gp["beta_L"], beta_S = gp["beta_S"],
                            psi_L = gp["psi_L"], psi_S = gp["psi_S"],
                            psi_LS = gp["psi_LS"])

  # ── Smriti: default (Pearson target, λ = 1.0) ────────────────────────────
  # NOTE: smriti reuses missForest's output as initial_imputation so the
  # Lagrangian routing benefit is measured in isolation.  time_sec captures
  # only the routing step; pipeline_time includes the full end-to-end cost
  # (missForest initialisation + smriti routing).
  time_sd <- system.time({
    imp_sd <- tryCatch(smriti_impute(df_miss, time_cols = 1:t_points,
                       initial_imputation = imp_mf, lambda = 1.0, robust = FALSE),
                       error = function(e) NULL)
    s_var_sd <- NA; s_se_sd <- NA; d_sd <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    if (!is.null(imp_sd)) {
      d_sd <- frob_dist(stats::cov(imp_sd[, 1:t_points]), true_cov)
      sv <- extract_slope_var(imp_sd, gcm_mod)
      s_var_sd <- sv["s_var"]; s_se_sd <- sv["s_se"]
      fit_sd <- tryCatch(growth(gcm_mod, data = imp_sd), error = function(e) NULL)
      if (!is.null(fit_sd)) gp <- extract_gcm_params(fit_sd)
    }
  })["elapsed"]
  res_list[[5]] <- make_row("Smriti_Default", d_sd, s_var_sd, s_se_sd,
                            unname(time_sd),
                            beta_L = gp["beta_L"], beta_S = gp["beta_S"],
                            psi_L = gp["psi_L"], psi_S = gp["psi_S"],
                            psi_LS = gp["psi_LS"],
                            pipeline_time = unname(time_mf) + unname(time_sd))

  # ── Smriti: robust (Spearman + MAD target, λ = 1.0) ──────────────────────
  time_sr <- system.time({
    imp_sr <- tryCatch(smriti_impute(df_miss, time_cols = 1:t_points,
                       initial_imputation = imp_mf, lambda = 1.0, robust = TRUE),
                       error = function(e) NULL)
    s_var_sr <- NA; s_se_sr <- NA; d_sr <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    if (!is.null(imp_sr)) {
      d_sr <- frob_dist(stats::cov(imp_sr[, 1:t_points]), true_cov)
      sv <- extract_slope_var(imp_sr, gcm_mod)
      s_var_sr <- sv["s_var"]; s_se_sr <- sv["s_se"]
      fit_sr <- tryCatch(growth(gcm_mod, data = imp_sr), error = function(e) NULL)
      if (!is.null(fit_sr)) gp <- extract_gcm_params(fit_sr)
    }
  })["elapsed"]
  res_list[[6]] <- make_row("Smriti_Robust", d_sr, s_var_sr, s_se_sr,
                            unname(time_sr),
                            beta_L = gp["beta_L"], beta_S = gp["beta_S"],
                            psi_L = gp["psi_L"], psi_S = gp["psi_S"],
                            psi_LS = gp["psi_LS"],
                            pipeline_time = unname(time_mf) + unname(time_sr))

  rm(df_true, df_miss, imp_mice, imp_mf, imp_mr, imp_sd, imp_sr)
  do.call(rbind, res_list)
}

# ── SLURM Array Dispatch ─────────────────────────────────────────────────────
conditions <- expand.grid(n = grid_n, miss = grid_miss, dist = grid_dist,
                          mech = grid_mech, stringsAsFactors = FALSE)
total_conditions <- nrow(conditions)

array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.na(array_id) && array_id >= 1 && array_id <= total_conditions) {
  current_conditions <- conditions[array_id, , drop = FALSE]
  output_file <- sprintf("sim_results/prod_results_%d.rds", array_id)
} else {
  current_conditions <- conditions
  output_file <- "sim_results/prod_results.rds"
}

# ── Execution ────────────────────────────────────────────────────────────────
cat(sprintf("Grid: %d conditions × %d reps × 6 methods = %d total rows\n",
            total_conditions, n_sims, total_conditions * n_sims * 6))
cat(sprintf("Parallel cores: %d\n", num_cores))
cat(sprintf("Output: %s\n\n", output_file))

dir.create("sim_results", showWarnings = FALSE)
dir.create("sim_raw_data", showWarnings = FALSE)

all_results <- list()

for (i in seq_len(nrow(current_conditions))) {
  params <- current_conditions[i, ]
  cond_idx <- if (!is.na(array_id)) array_id else i
  cat(sprintf("[%s] Condition %d/%d: N=%d, Miss=%.2f, Dist=%s, Mech=%s\n",
              Sys.time(), cond_idx, total_conditions, params$n, params$miss,
              params$dist, params$mech))

  out_list <- mclapply(seq_len(n_sims), function(s) run_iteration(s, params),
                        mc.cores = num_cores, mc.preschedule = TRUE)

  valid <- out_list[sapply(out_list, is.data.frame)]
  if (length(valid) > 0) {
    cond_df <- do.call(rbind, valid)
    all_results[[i]] <- cond_df
    saveRDS(do.call(rbind, all_results), output_file)
  }

  rm(out_list, valid)
  gc()
}

cat("\nProduction simulation complete. Results saved to ", output_file, "\n")
