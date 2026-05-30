library(MASS)
library(smriti)
library(parallel)
library(lavaan)
library(missForest)

# Neutralize multi-threaded BLAS to prevent resource contention during forking
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

set.seed(20250523)

# ── Speed Toggle (must precede grid + core count) ────────────────────────────
# "coarse"  — quick pass (~1 hr local):   3 N, 5 lambda, 2 robust, 2 dist, 100 reps
# "full"    — proper tuning (~6 hr local): 5 N, 7 lambda, 2 robust, 2 dist, 200 reps
tune_mode <- "coarse"
if (tune_mode == "coarse") {
  grid_n      <- c(200, 500, 5000)
  grid_lambda <- c(0.01, 0.05, 0.1, 1.0, 5.0)
  n_sims      <- 100
  grid_dist   <- c("Lognormal", "Normal")
} else {
  grid_n      <- c(100, 200, 500, 1000, 5000)
  grid_lambda <- c(0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0)
  n_sims      <- 200
  grid_dist   <- c("Lognormal", "Normal")
}

# ── Core Allocation ──────────────────────────────────────────────────────────
# mclapply forks lightweight R processes; peak RAM per worker at N=5000
# is ~300 MB, so 14 workers on 16 GB is safe.
num_cores <- if (tune_mode == "coarse") {
  min(14, parallel::detectCores() - 1)
} else {
  min(12, parallel::detectCores() - 1)
}

# ── Tuning Grid ──────────────────────────────────────────────────────────────
# FIML is the sole baseline — MICE, missForest, missRanger are excluded.
# The smriti Lagrangian step is ~0.02s per variant, so sweeping 5 lambda ×
# 2 robust combinations adds only ~15% overhead over a single config.
# grid_n, grid_lambda, n_sims, and grid_dist are set by the speed toggle above.
grid_miss   <- c(0.05, 0.15, 0.30)
grid_mech   <- c("MAR", "MNAR")
t_points    <- 4

# Smriti hyperparameter candidates to evaluate against FIML
grid_robust <- c(TRUE, FALSE)

# ── True Population Parameters ───────────────────────────────────────────────
mu_i <- 6.0; mu_s <- 2.0; v_i <- 1.0; v_s <- 1.0; c_is <- 0.0; v_e <- 1.0

# Compute True Population Covariance Matrix
true_cov <- matrix(0, t_points, t_points)
for (r in 1:t_points) {
  for (c in 1:t_points) {
    true_cov[r, c] <- v_i + (r - 1) * (c - 1) * v_s + ((r - 1) + (c - 1)) * c_is
    if (r == c) true_cov[r, c] <- true_cov[r, c] + v_e
  }
}

# ── Metrics ──────────────────────────────────────────────────────────────────
# Primary tuning objective: Frobenius distance from true covariance
frob_dist <- function(m1, m2) sqrt(sum((m1 - m2)^2))

# Secondary diagnostic: relative bias (%) of a scalar estimate
rel_bias <- function(est, truth) 100 * (est - truth) / truth

# ── Data Generation Engine ───────────────────────────────────────────────────
generate_data <- function(n, dist) {
  latent_vars <- mvrnorm(n, mu = c(mu_i, mu_s),
                         Sigma = matrix(c(v_i, c_is, c_is, v_s), 2, 2))

  if (dist == "Lognormal") {
    # Induce lognormality while preserving target moments
    latent_vars <- exp(latent_vars)
    latent_vars[,1] <- scale(latent_vars[,1]) * sqrt(v_i) + mu_i
    latent_vars[,2] <- scale(latent_vars[,2]) * sqrt(v_s) + mu_s
  }

  data_mat <- matrix(0, n, t_points)
  for (j in 1:t_points) {
    # Error term: log-normal for the Lognormal condition, Gaussian otherwise
    err <- if (dist == "Lognormal") {
      scale(exp(rnorm(n))) * sqrt(v_e)
    } else {
      rnorm(n, 0, sqrt(v_e))
    }
    data_mat[, j] <- latent_vars[, 1] + (j - 1) * latent_vars[, 2] + err
  }

  if (dist == "Outlier") {
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
    # Probabilistic Logistic MAR (safe against propagating NAs)
    for (t in 1:(t_points - 1)) {
      idx <- which(!is.na(df_miss[, t]))
      if (length(idx) > 0) {
        x_prev <- scale(df_miss[idx, t])
        p_miss <- 1 / (1 + exp(-(x_prev - qnorm(1 - rate))))
        drop_idx <- idx[rbinom(length(idx), 1, p_miss) == 1]
        df_miss[drop_idx, (t + 1)] <- NA
      }
    }
  } else if (mech == "MNAR") {
    # Latent-dependent dropout
    cor_ab <- 0.8; a <- cor_ab / sqrt(1 - cor_ab^2)
    aux_var <- a * df$true_slope + rnorm(n, 0, 1)
    miss_rate_t <- 2 * rate / (t_points - 1)
    for (j in 2:t_points) {
      crit <- qnorm((1 - (j - 1) * miss_rate_t), mean = a * mu_s, sd = sqrt(a^2 + 1))
      df_miss[which(aux_var > crit), j] <- NA
    }
  }

  df_miss$true_slope <- NULL
  return(df_miss)
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
  res_list <- list()

  df_true <- generate_data(params$n, params$dist)
  df_miss <- apply_missingness(df_true, params$miss, params$mech)

  # ── FIML Baseline (run once per replication) ─────────────────────────────
  time_fiml <- system.time({
    fit_fiml <- tryCatch(
      growth(gcm_mod, data = df_miss, missing = "fiml"),
      error = function(e) NULL
    )
    s_var_f <- NA; s_se_f <- NA; d_fiml <- NA; rb_f <- NA
    if (!is.null(fit_fiml)) {
      pt <- parameterEstimates(fit_fiml)
      row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
      if (nrow(row) > 0) {
        s_var_f <- row$est[1]
        s_se_f  <- row$se[1]
        rb_f    <- rel_bias(s_var_f, v_s)
      }
      # Model-implied covariance to avoid lavPredict conditional-expectation bias
      implied_cov <- tryCatch(
        lavaan::lavInspect(fit_fiml, "cov.ov"),
        error = function(e) NULL
      )
      if (!is.null(implied_cov)) d_fiml <- frob_dist(implied_cov, true_cov)
    }
  })["elapsed"]
  res_list[[1]] <- data.frame(
    sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist,
    mech = params$mech, method = "FIML", lambda = NA, robust = NA,
    f_dist = d_fiml, s_var = s_var_f, s_se = s_se_f, rel_bias = rb_f,
    time_sec = unname(time_fiml)
  )

  # ── Initial Imputation (shared across smriti variants) ──────────────────
  imp_mf <- tryCatch(
    missForest::missForest(df_miss, verbose = FALSE)$ximp,
    error = function(e) NULL
  )
  if (is.null(imp_mf)) {
    # Fallback: column-mean initialisation
    imp_mf <- as.matrix(df_miss[, 1:t_points])
    for (j in 1:t_points) {
      na_idx <- is.na(imp_mf[, j])
      if (any(na_idx)) imp_mf[na_idx, j] <- mean(imp_mf[, j], na.rm = TRUE)
    }
  }

  # ── Smriti variants (grid over lambda × robust) ─────────────────────────
  lambdas <- params$lambdas[[1]]
  robusts <- params$robusts[[1]]
  for (lam in lambdas) {
    for (rb in robusts) {
      tag <- sprintf("Smriti_l%.2f_%s", lam, if (rb) "R" else "S")
      time_sm <- system.time({
        imp_sm <- tryCatch(
          smriti_impute(df_miss, time_cols = 1:t_points,
                        initial_imputation = imp_mf,
                        lambda = lam, robust = rb),
          error = function(e) NULL
        )
        s_var_s <- NA; s_se_s <- NA; d_sm <- NA; rb_s <- NA
        if (!is.null(imp_sm)) {
          cov_imp <- stats::cov(imp_sm[, 1:t_points])
          d_sm <- frob_dist(cov_imp, true_cov)
          fit_sm <- tryCatch(
            growth(gcm_mod, data = imp_sm),
            error = function(e) NULL
          )
          if (!is.null(fit_sm)) {
            pt <- parameterEstimates(fit_sm)
            row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
            if (nrow(row) > 0) {
              s_var_s <- row$est[1]
              s_se_s  <- row$se[1]
              rb_s    <- rel_bias(s_var_s, v_s)
            }
          }
        }
      })["elapsed"]
      res_list[[length(res_list) + 1]] <- data.frame(
        sim_id = sim_id, N = params$n, miss = params$miss, dist = params$dist,
        mech = params$mech, method = tag, lambda = lam, robust = rb,
        f_dist = d_sm, s_var = s_var_s, s_se = s_se_s, rel_bias = rb_s,
        time_sec = unname(time_sm)
      )
    }
  }

  rm(df_true, df_miss, imp_mf)
  do.call(rbind, res_list)
}

# ── Condition / Hyperparameter Expansion ─────────────────────────────────────
conditions <- expand.grid(
  n = grid_n, miss = grid_miss, dist = grid_dist, mech = grid_mech,
  stringsAsFactors = FALSE
)
# Each condition carries the full hyperparameter grid for smriti
conditions$lambdas <- list(grid_lambda)
conditions$robusts <- list(grid_robust)
total_conditions <- nrow(conditions)
n_variants <- 1 + length(grid_lambda) * length(grid_robust)  # FIML + smriti combos

# ── SLURM Array Dispatch ─────────────────────────────────────────────────────
# When running under a SLURM job array, each task processes exactly one
# condition.  Locally the full grid is processed sequentially.
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.na(array_id) && array_id >= 1 && array_id <= total_conditions) {
  current_conditions <- conditions[array_id, , drop = FALSE]
  output_file <- sprintf("sim_results/tune_results_%d.rds", array_id)
} else {
  current_conditions <- conditions
  output_file <- "sim_results/tune_results.rds"
}

# ── Execution ────────────────────────────────────────────────────────────────
cat(sprintf("Speed Toggle: %s\n", tune_mode))
cat(sprintf("Grid: %d conditions × %d reps × %d variants = %d total rows\n",
            total_conditions, n_sims, n_variants,
            total_conditions * n_sims * n_variants))
cat(sprintf("Lambda grid: %s\n", paste(grid_lambda, collapse = ", ")))
cat(sprintf("Parallel cores: %d\n", num_cores))
cat(sprintf("Output: %s\n\n", output_file))

dir.create("sim_results", showWarnings = FALSE)

all_results <- list()

for (i in seq_len(nrow(current_conditions))) {
  params <- current_conditions[i, ]
  cond_idx <- if (!is.na(array_id)) array_id else i
  cat(sprintf("[%s] Condition %d/%d: N=%d, Miss=%.2f, Dist=%s, Mech=%s\n",
              Sys.time(), cond_idx, total_conditions, params$n, params$miss,
              params$dist, params$mech))

  out_list <- mclapply(
    seq_len(n_sims),
    function(s) run_iteration(s, params),
    mc.cores = num_cores,
    mc.preschedule = TRUE
  )

  valid <- out_list[sapply(out_list, is.data.frame)]
  if (length(valid) > 0) {
    cond_df <- do.call(rbind, valid)
    all_results[[i]] <- cond_df
    saveRDS(do.call(rbind, all_results), output_file)
  }

  rm(out_list, valid)
  gc()
}

cat(sprintf("\nTuning complete. Results saved to %s\n", output_file))
