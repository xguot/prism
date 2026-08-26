# Rivanna local library path (R_LIBS_USER may be ignored by site .Renviron)
.libPaths(c("~/R/rivanna-lib", .libPaths()))

library(MASS)
library(prism)
library(parallel)
library(lavaan)
library(mice)
library(missForest)

# Neutralize multi-threaded BLAS to prevent resource contention
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# Core allocation for production run
num_cores <- min(24, parallel::detectCores() - 1)
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
seed_base <- 20250523
set.seed(if (is.na(array_id)) seed_base else seed_base + array_id)

# ── Simulation Grid ───────────────────────────────────────────────────────────
# Aligned with Tang & Tong (UVA) manuscript:
#   model   — GCM (latent growth) and SEM (two-factor structural regression)
#   N       — 100, 200, 500, 1000, 5000, 10000
#   miss    — 5%, 10%, 15%, 30%
#   dist    — Normal, t(5), Outlier, Lognormal
#   mech    — MAR, MNAR
grid_model <- c("GCM", "SEM")
grid_n    <- c(100, 200, 500, 1000, 5000, 10000)
grid_miss <- c(0.05, 0.10, 0.15, 0.30)
grid_dist <- c("Normal", "t5", "Outlier", "Lognormal")
grid_mech <- c("MAR", "MNAR")
n_sims    <- 500
t_points  <- 4

# Two-level PRISM_MI: when TRUE, level-1 missing-value uncertainty comes from
# the default bootstrap-weighted forest initializer (proper two-level MI);
# when FALSE, the deterministic missForest matrix is reused for every draw
# (v1-style level-1; cheaper and directly comparable to the PRISM row).
two_level_mi <- FALSE
options(prism.num_trees = 100L)   # tractable forest size for the study

# ── GCM: True Population Parameters ───────────────────────────────────────────
#   σ²_i = 1, σ²_s = 1, σ_is = 0, σ²_e = 1, μ_i = 6, μ_s = 2
mu_i <- 6.0; mu_s <- 2.0; v_i <- 1.0; v_s <- 1.0; c_is <- 0.0; v_e <- 1.0

gcm_mod <- "
  i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
  s =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
  i ~~ s
  i ~~ i
  s ~~ s
"

# True population covariance matrix (4 × 4)
true_cov <- matrix(0, t_points, t_points)
for (r in 1:t_points) {
  for (c in 1:t_points) {
    true_cov[r, c] <- v_i + (r - 1) * (c - 1) * v_s + ((r - 1) + (c - 1)) * c_is
    if (r == c) true_cov[r, c] <- true_cov[r, c] + v_e
  }
}

# ── SEM: True Population Parameters ───────────────────────────────────────────
#   Two-factor model with a structural regression eta2 ~ eta1:
#     loadings λ = (0.7, 0.75, 0.8) per factor, b21 = 0.5,
#     factor variances ψ_F1 = ψ_F2 = 1, indicator variances 1 - λ².
sem_mod <- "
  eta1 =~ y1 + y2 + y3
  eta2 =~ y4 + y5 + y6
  eta2 ~ eta1
"

sem_lam     <- c(0.7, 0.75, 0.8)
b21_true    <- 0.5
psi_F1_true <- 1.0
# psi_F2 is the RESIDUAL variance of the endogenous factor:
# eta2 = b21 * eta1 + zeta with Var(zeta) = 1 - b21^2 * psi_F1
psi_F2_true <- 1.0 - b21_true^2 * psi_F1_true

true_cov_sem <- (function() {
  Sigma <- matrix(0, 6, 6)
  Sigma[1:3, 1:3] <- outer(sem_lam, sem_lam) + diag(1 - sem_lam^2)
  Sigma[4:6, 4:6] <- outer(sem_lam, sem_lam) + diag(1 - sem_lam^2)
  cross <- b21_true * outer(sem_lam, sem_lam)
  Sigma[4:6, 1:3] <- cross
  Sigma[1:3, 4:6] <- t(cross)
  Sigma
})()

# ── Metrics ──────────────────────────────────────────────────────────────────
frob_dist <- function(m1, m2) sqrt(sum((m1 - m2)^2))
rel_bias  <- function(est, truth) {
  if (abs(truth) < 1e-12) return(est - truth)  # Tang & Tong: RB = raw bias when θ = 0
  100 * (est - truth) / truth
}

# ── Diagnostic helpers (mean drift and PRISM diagnostics) ─────────────────────
# delta_mu: projection drift — max |mean(completed) - mean(initial X0)|
# mean_bias: max |mean(completed) - mean(true complete data)|
delta_mu_fn  <- function(comp, x0) max(abs(colMeans(comp) - colMeans(x0)))
mean_bias_fn <- function(comp, truth) max(abs(colMeans(comp) - colMeans(truth)))

mi_mean <- function(fn, imp_list, ref) {
  mean(vapply(imp_list, function(x) fn(x, ref), numeric(1)))
}

prism_diag_fn <- function(imp) {
  d <- attr(imp, "prism_diagnostics")
  if (is.null(d)) {
    return(list(r_kkT = NA_real_, feas_gap = NA_real_, fidelity = NA_real_,
                status = NA_character_))
  }
  list(r_kkT    = d$r_kkT,
       feas_gap = d$feas_gap,
       fidelity = d$fidelity,
       status   = d$status)
}

mi_prism_diag_fn <- function(imp_list) {
  ds <- lapply(imp_list, prism_diag_fn)
  list(
    r_kkT    = mean(vapply(ds, `[[`, numeric(1), "r_kkT")),
    feas_gap = mean(vapply(ds, `[[`, numeric(1), "feas_gap")),
    fidelity = mean(vapply(ds, `[[`, numeric(1), "fidelity")),
    status   = if (length(ds) > 0) ds[[1]]$status else NA_character_
  )
}

# ── Data Generation Engines ───────────────────────────────────────────────────
generate_gcm_data <- function(n, dist) {
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

    if (dist == "Outlier") {
      # Occasion-specific measurement error outliers (5%)
      out_idx <- sample(seq_len(n), floor(0.05 * n))
      err[out_idx] <- err[out_idx] + 5.0
    }

    data_mat[, j] <- latent_vars[, 1] + (j - 1) * latent_vars[, 2] + err
  }

  df <- as.data.frame(data_mat)
  colnames(df) <- paste0("T", 1:t_points)
  df$true_slope <- latent_vars[, 2]
  return(df)
}

generate_sem_data <- function(n, dist) {
  eta1 <- rnorm(n, 0, 1)
  eta2 <- b21_true * eta1 + rnorm(n, 0, sqrt(1 - b21_true^2))

  gen_block <- function(eta) {
    vapply(seq_along(sem_lam), function(k) {
      err <- if (dist == "Lognormal") {
        scale(exp(rnorm(n))) * sqrt(1 - sem_lam[k]^2)
      } else if (dist == "t5") {
        rt(n, df = 5) * sqrt((1 - sem_lam[k]^2) * (5 - 2) / 5)
      } else {
        rnorm(n, 0, sqrt(1 - sem_lam[k]^2))
      }

      if (dist == "Outlier") {
        out_idx <- sample(seq_len(n), floor(0.05 * n))
        err[out_idx] <- err[out_idx] + 5.0
      }

      sem_lam[k] * eta + err
    }, numeric(n))
  }

  df <- data.frame(gen_block(eta1), gen_block(eta2))
  colnames(df) <- paste0("y", 1:6)
  df
}

# ── Missingness Engines ───────────────────────────────────────────────────────
# GCM convention (Tang & Tong 2025): dropout at T_{t+1} driven by low T_t
# under MAR; latent-slope-driven dropout under MNAR.
apply_missingness <- function(df, rate, mech) {
  df_miss <- df; n <- nrow(df)

  if (mech == "MAR") {
    # Deterministic threshold MAR (Tang & Tong 2025 convention)
    # Low values on T_t cause missingness on T_{t+1} (monotone dropout).
    # Drop budget accumulates linearly: t * miss_n_per_step rows at step t.
    # With T = 4 the realized cellwise rate is 10 * miss_n_per_step / (4n),
    # i.e. about (5/3) * rate; actual_miss is recorded per replication.
    miss_n_per_step <- round(2 * n * rate / (t_points - 1))

    for (t in 1:(t_points - 1)) {
      # Select drop targets only among subjects still observed at time t:
      # subjects dropped at earlier steps already carry NA and must not
      # consume this step's drop budget
      obs_idx <- which(!is.na(df_miss[, t]))
      if (length(obs_idx) > 0) {
        # smallest values of T_t first, restricted to the still-observed rows
        order_idx <- obs_idx[order(df_miss[obs_idx, t], decreasing = TRUE)]
        n_take    <- min(t * miss_n_per_step, length(obs_idx))
        drop_idx  <- tail(order_idx, n_take)
        df_miss[drop_idx, (t + 1):t_points] <- NA # Dropout: once missing, stay missing
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

# SEM conventions: `rate` is the per-variable missingness fraction.
#   MAR  — chain mechanism: low values of y_{j-1} drive missingness in y_j
#   MNAR — self-missingness: the lowest values of y_j are missing in y_j
apply_missingness_sem <- function(df, rate, mech) {
  df_miss <- df; n <- nrow(df); p <- ncol(df)

  if (mech == "MAR") {
    for (j in 2:p) {
      obs_idx <- which(!is.na(df_miss[, j - 1]))
      if (length(obs_idx) == 0) next
      order_idx <- obs_idx[order(df_miss[obs_idx, j - 1], decreasing = TRUE)]
      n_miss_j  <- round(length(obs_idx) * rate)
      drop_idx  <- tail(order_idx, n_miss_j)
      df_miss[drop_idx, j] <- NA
    }
  } else if (mech == "MNAR") {
    for (j in 1:p) {
      obs_idx <- which(!is.na(df_miss[, j]))
      order_idx <- obs_idx[order(df_miss[obs_idx, j], decreasing = TRUE)]
      n_miss_j  <- round(length(obs_idx) * rate)
      drop_idx  <- tail(order_idx, n_miss_j)
      df_miss[drop_idx, j] <- NA
    }
  }

  list(data = df_miss, actual_miss = mean(is.na(df_miss)))
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

# ── Helper: extract SEM parameters from lavaan fit ────────────────────────────
# Returns a named vector of 3 structural parameters:
#   b21 (structural path), psi_F1, psi_F2 (factor variances)
extract_sem_params <- function(fit, type = "est") {
  pt <- parameterEstimates(fit)
  get_val <- function(lhs, op, rhs) {
    row <- pt[pt$lhs == lhs & pt$op == op & pt$rhs == rhs, ]
    if (nrow(row) > 0) row[[type]][1] else NA_real_
  }
  c(
    b21    = get_val("eta2", "~", "eta1"),
    psi_F1 = get_val("eta1", "~~", "eta1"),
    psi_F2 = get_val("eta2", "~~", "eta2")
  )
}

# ── Helper: extract structural path + SE (SEM key metric) ─────────────────────
# Null fits (failed optimization) yield NA instead of an error, so a single
# failed refit cannot abort the whole replication.
extract_sem_path <- function(fit) {
  if (is.null(fit)) return(c(b = NA, b_se = NA))
  pt <- parameterEstimates(fit)
  row <- pt[pt$lhs == "eta2" & pt$op == "~" & pt$rhs == "eta1", ]
  if (nrow(row) > 0) c(b = row$est[1], b_se = row$se[1])
  else c(b = NA, b_se = NA)
}

# ── Shared result schema ──────────────────────────────────────────────────────
# One row per (simulation, model, method).  GCM-specific parameter columns are
# NA for SEM rows and vice versa; prism diagnostics and mean-drift columns are
# NA for methods that do not produce them.
make_result_row <- function(sim_id, params, act_m, model, method, f_dist,
                            time_sec, pipeline_time = time_sec,
                            converged = NA_integer_,
                            delta_mu = NA_real_, mean_bias = NA_real_,
                            r_kkT = NA_real_, feas_gap = NA_real_,
                            fidelity = NA_real_, status = NA_character_,
                            s_var = NA_real_, s_var_bias = NA_real_,
                            s_se = NA_real_,
                            est_L = NA_real_, est_S = NA_real_,
                            est_var_L = NA_real_, est_var_S = NA_real_,
                            est_cov_LS = NA_real_,
                            se_L = NA_real_, se_S = NA_real_,
                            se_var_L = NA_real_, se_var_S = NA_real_,
                            se_cov_LS = NA_real_,
                            bias_L = NA_real_, bias_S = NA_real_,
                            bias_var_L = NA_real_, bias_var_S = NA_real_,
                            bias_cov_LS = NA_real_,
                            est_b21 = NA_real_, se_b21 = NA_real_,
                            bias_b21 = NA_real_,
                            est_var_F1 = NA_real_, se_var_F1 = NA_real_,
                            bias_var_F1 = NA_real_,
                            est_var_F2 = NA_real_, se_var_F2 = NA_real_,
                            bias_var_F2 = NA_real_) {
  # diagnostics list fields may be absent (NULL) for methods without them
  nz <- function(x) if (is.null(x)) NA else x
  data.frame(
    sim_id = sim_id, model = model, N = params$n, miss = params$miss,
    actual_miss = act_m,
    dist = params$dist, mech = params$mech,
    method       = method,
    f_dist       = f_dist,
    converged    = converged,
    delta_mu     = nz(delta_mu),
    mean_bias    = nz(mean_bias),
    r_kkT        = nz(r_kkT),
    feas_gap     = nz(feas_gap),
    fidelity     = nz(fidelity),
    status       = nz(status),
    s_var        = s_var,
    s_var_bias   = s_var_bias,
    s_se         = s_se,
    est_L        = est_L,       est_S = est_S,
    est_var_L    = est_var_L,   est_var_S = est_var_S,
    est_cov_LS   = est_cov_LS,
    se_L         = se_L,        se_S = se_S,
    se_var_L     = se_var_L,    se_var_S = se_var_S,
    se_cov_LS    = se_cov_LS,
    bias_L       = bias_L,      bias_S = bias_S,
    bias_var_L   = bias_var_L,  bias_var_S = bias_var_S,
    bias_cov_LS  = bias_cov_LS,
    est_b21      = est_b21,     se_b21 = se_b21,
    bias_b21     = bias_b21,
    est_var_F1   = est_var_F1,  se_var_F1 = se_var_F1,
    bias_var_F1  = bias_var_F1,
    est_var_F2   = est_var_F2,  se_var_F2 = se_var_F2,
    bias_var_F2  = bias_var_F2,
    time_sec     = time_sec,
    pipeline_time = pipeline_time,
    stringsAsFactors = FALSE
  )
}

# ── GCM Single-Replication Worker ─────────────────────────────────────────────
run_iteration_gcm <- function(sim_id, params) {
  mk <- function(method, f_dist, s_var, s_se, time_sec, gp, gs,
                 pipeline_time = time_sec, diag = list()) {
    make_result_row(
      sim_id = sim_id, params = params, act_m = act_m, model = "GCM",
      method = method, f_dist = f_dist,
      time_sec = unname(time_sec), pipeline_time = unname(pipeline_time),
      converged = as.integer(!is.na(s_var) & is.finite(s_var) &&
                             !is.na(gp["beta_L"]) & is.finite(gp["beta_L"])),
      s_var = s_var, s_var_bias = rel_bias(s_var, v_s), s_se = s_se,
      est_L = gp["beta_L"], est_S = gp["beta_S"],
      est_var_L = gp["psi_L"], est_var_S = gp["psi_S"],
      est_cov_LS = gp["psi_LS"],
      se_L = gs["beta_L"], se_S = gs["beta_S"],
      se_var_L = gs["psi_L"], se_var_S = gs["psi_S"],
      se_cov_LS = gs["psi_LS"],
      bias_L = rel_bias(gp["beta_L"], mu_i),
      bias_S = rel_bias(gp["beta_S"], mu_s),
      bias_var_L = rel_bias(gp["psi_L"], v_i),
      bias_var_S = rel_bias(gp["psi_S"], v_s),
      bias_cov_LS = rel_bias(gp["psi_LS"], c_is),
      delta_mu = diag$delta_mu, mean_bias = diag$mean_bias,
      r_kkT = diag$r_kkT, feas_gap = diag$feas_gap,
      fidelity = diag$fidelity, status = diag$status
    )
  }
  res_list <- list()

  df_true <- generate_gcm_data(params$n, params$dist)
  miss_out <- apply_missingness(df_true, params$miss, params$mech)
  df_miss <- miss_out$data
  act_m   <- miss_out$actual_miss
  df_true_ov <- df_true[, 1:t_points]

  # ── FIML Baseline ────────────────────────────────────────────────────────
  time_fiml <- system.time({
    fit_fiml <- tryCatch(growth(gcm_mod, data = df_miss, missing = "fiml"),
                         error = function(e) NULL)
    s_var_f <- NA; s_se_f <- NA; d_fiml <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    gs <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    if (!is.null(fit_fiml)) {
      pt <- parameterEstimates(fit_fiml)
      row <- pt[pt$lhs == "s" & pt$op == "~~" & pt$rhs == "s", ]
      if (nrow(row) > 0) { s_var_f <- row$est[1]; s_se_f <- row$se[1] }
      implied_cov <- tryCatch(lavaan::lavInspect(fit_fiml, "cov.ov"),
                              error = function(e) NULL)
      if (!is.null(implied_cov)) d_fiml <- frob_dist(implied_cov, true_cov)
      gp <- extract_gcm_params(fit_fiml)
      gs <- extract_gcm_params(fit_fiml, "se")
    }
  })["elapsed"]
  res_list[[1]] <- mk("FIML", d_fiml, s_var_f, s_se_f, time_fiml, gp, gs)

  # ── FIML + lavPredict: conditional-expectation completed data ─────────────
  # lavPredict(type="ov") gives the conditional expectation of each missing
  # value given the observed data and the FIML parameter estimates.  This is
  # the "obvious" way to get completed data from a FIML model — but the
  # conditional expectations are shrunk toward the mean, attenuating variance.
  # Included here as a baseline to demonstrate that the naive FIML-completed
  # dataset fails at the covariance level, motivating the prism projection.
  time_fp <- system.time({
    s_var_fp <- NA; s_se_fp <- NA; d_fp <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    gs <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    fp_diag <- list()
    if (!is.null(fit_fiml)) {
      imp_fp <- tryCatch(lavaan::lavPredict(fit_fiml, type = "ov"),
                         error = function(e) NULL)
      if (!is.null(imp_fp)) {
        d_fp <- frob_dist(stats::cov(imp_fp[, 1:t_points]), true_cov)
        colnames(imp_fp) <- paste0("T", 1:t_points)
        sv <- extract_slope_var(as.data.frame(imp_fp), gcm_mod)
        s_var_fp <- sv["s_var"]; s_se_fp <- sv["s_se"]
        fit_fp <- tryCatch(growth(gcm_mod, data = as.data.frame(imp_fp)),
                           error = function(e) NULL)
        if (!is.null(fit_fp)) { gp <- extract_gcm_params(fit_fp); gs <- extract_gcm_params(fit_fp, "se") }
        fp_diag <- list(mean_bias = mean_bias_fn(imp_fp, df_true_ov))
      }
    }
  })["elapsed"]
  res_list[[2]] <- mk("FIML_lavPredict", d_fp, s_var_fp, s_se_fp, time_fp,
                      gp, gs, diag = fp_diag)

  # ── MICE (random forest via ranger, MI m=20) ──────────────────────────────
  time_mice_rf <- system.time({
    m_mice_rf <- 20
    imp_mice_rf_list <- tryCatch({
      imp_obj_rf <- mice::mice(df_miss, m = m_mice_rf, method = "rf", printFlag = FALSE)
      mice::complete(imp_obj_rf, "all")
    }, error = function(e) NULL)

    s_var_mr <- NA; s_se_mr <- NA; d_mr <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    gps <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    mr_diag <- list()

    if (!is.null(imp_mice_rf_list)) {
      cov_list <- lapply(imp_mice_rf_list, function(x) stats::cov(x[, 1:t_points]))
      avg_cov <- Reduce("+", cov_list) / length(cov_list)
      d_mr <- frob_dist(avg_cov, true_cov)
      mr_diag <- list(mean_bias = mi_mean(mean_bias_fn, imp_mice_rf_list, df_true_ov))

      fit_results <- lapply(imp_mice_rf_list, function(ds) {
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
        s_var_mr <- gp["psi_S"]; s_se_mr <- gps["psi_S"]
      }
    }
  })["elapsed"]
  res_list[[3]] <- mk("MICE", d_mr, s_var_mr, s_se_mr, time_mice_rf,
                      gp, gps, diag = mr_diag)

  # ── missForest Baseline ───────────────────────────────────────────────────
  time_mf <- system.time({
    imp_mf <- tryCatch(missForest::missForest(df_miss, verbose = FALSE)$ximp,
                       error = function(e) NULL)
    s_var_mf <- NA; s_se_mf <- NA; d_mf <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    gs <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    mf_diag <- list()
    if (!is.null(imp_mf)) {
      d_mf <- frob_dist(stats::cov(imp_mf[, 1:t_points]), true_cov)
      sv <- extract_slope_var(imp_mf, gcm_mod)
      s_var_mf <- sv["s_var"]; s_se_mf <- sv["s_se"]
      fit_mf <- tryCatch(growth(gcm_mod, data = imp_mf), error = function(e) NULL)
      if (!is.null(fit_mf)) { gp <- extract_gcm_params(fit_mf); gs <- extract_gcm_params(fit_mf, "se") }
      mf_diag <- list(mean_bias = mean_bias_fn(imp_mf, df_true_ov))
    }
  })["elapsed"]
  res_list[[4]] <- mk("missForest", d_mf, s_var_mf, s_se_mf, time_mf,
                      gp, gs, diag = mf_diag)

  # ── prism_fiml: FIML model-implied Σ target ───────────────────────────
  # Uses prism_fiml() which fits a lavaan growth model with FIML to extract
  # the model-implied covariance as the structural target, then projects the
  # missForest initial imputation toward it.
  time_sf <- system.time({
    imp_sf <- tryCatch(suppressWarnings(
      prism_fiml(df_miss, model = gcm_mod,
                 initial_imputation = imp_mf, lambda_sigma = 1.0)),
      error = function(e) NULL)
    s_var_sf <- NA; s_se_sf <- NA; d_sf <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    gs <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    sf_diag <- list()
    if (!is.null(imp_sf)) {
      d_sf <- frob_dist(stats::cov(imp_sf[, 1:t_points]), true_cov)
      sv <- extract_slope_var(imp_sf, gcm_mod)
      s_var_sf <- sv["s_var"]; s_se_sf <- sv["s_se"]
      fit_sf <- tryCatch(growth(gcm_mod, data = imp_sf), error = function(e) NULL)
      if (!is.null(fit_sf)) { gp <- extract_gcm_params(fit_sf); gs <- extract_gcm_params(fit_sf, "se") }
      sf_diag <- c(list(delta_mu = if (!is.null(imp_mf)) delta_mu_fn(imp_sf[, 1:t_points], imp_mf[, 1:t_points]) else NA_real_,
                        mean_bias = mean_bias_fn(imp_sf[, 1:t_points], df_true_ov)),
                   prism_diag_fn(imp_sf))
    }
  })["elapsed"]
  res_list[[5]] <- mk("PRISM", d_sf, s_var_sf, s_se_sf, time_sf, gp, gs,
                      pipeline_time = unname(time_mf) + unname(time_sf),
                      diag = sf_diag)

  # ── prism_mi: Proper Multiple Imputation ────────────────────────────────────
  time_smi <- system.time({
    m_prism <- 20
    imp_smi_list <- tryCatch(suppressWarnings({
      fit_fiml_base <- growth(gcm_mod, data = df_miss, missing = "fiml")
      if (two_level_mi) {
        prism_mi(df_miss, fit_fiml_base, m = m_prism, lambda_sigma = 1.0)
      } else {
        prism_mi(df_miss, fit_fiml_base, m = m_prism,
                 initial_imputation = imp_mf, lambda_sigma = 1.0)
      }
    }), error = function(e) NULL)

    s_var_smi <- NA; s_se_smi <- NA; d_smi <- NA
    gp <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    gps <- c(beta_L = NA, beta_S = NA, psi_L = NA, psi_S = NA, psi_LS = NA)
    smi_diag <- list()

    if (!is.null(imp_smi_list) && length(imp_smi_list) == m_prism) {
      cov_list <- lapply(imp_smi_list, function(x) stats::cov(x[, 1:t_points]))
      avg_cov <- Reduce("+", cov_list) / length(cov_list)
      d_smi <- frob_dist(avg_cov, true_cov)
      smi_diag <- c(list(delta_mu  = mi_mean(delta_mu_fn,  imp_smi_list, imp_mf[, 1:t_points]),
                         mean_bias = mi_mean(mean_bias_fn, imp_smi_list, df_true_ov)),
                    mi_prism_diag_fn(imp_smi_list))

      fit_results <- lapply(imp_smi_list, function(ds) {
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
        s_var_smi <- gp["psi_S"]; s_se_smi <- gps["psi_S"]
      }
    }
  })["elapsed"]
  res_list[[6]] <- mk("PRISM_MI", d_smi, s_var_smi, s_se_smi, time_smi,
                      gp, gps,
                      pipeline_time = unname(time_mf) + unname(time_smi),
                      diag = smi_diag)

  rm(df_true, df_miss, imp_mice_rf_list, imp_mf, imp_sf, imp_fp, imp_smi_list)
  do.call(rbind, res_list)
}

# ── SEM Single-Replication Worker ─────────────────────────────────────────────
# Mirrors the GCM worker: the four baselines are model-agnostic, prism_sem /
# prism_mi replace prism_fiml / prism_mi, and the key parameter is the
# structural path b21 (analogous to the GCM slope variance).
run_iteration_sem <- function(sim_id, params) {
  mk <- function(method, f_dist, b_est, b_se, time_sec, gp, gs,
                 pipeline_time = time_sec, diag = list()) {
    make_result_row(
      sim_id = sim_id, params = params, act_m = act_m, model = "SEM",
      method = method, f_dist = f_dist,
      time_sec = unname(time_sec), pipeline_time = unname(pipeline_time),
      converged = as.integer(!is.na(b_est) & is.finite(b_est) &&
                             !is.na(gp["b21"]) & is.finite(gp["b21"])),
      est_b21 = gp["b21"], se_b21 = gs["b21"],
      bias_b21 = rel_bias(gp["b21"], b21_true),
      est_var_F1 = gp["psi_F1"], se_var_F1 = gs["psi_F1"],
      bias_var_F1 = rel_bias(gp["psi_F1"], psi_F1_true),
      est_var_F2 = gp["psi_F2"], se_var_F2 = gs["psi_F2"],
      bias_var_F2 = rel_bias(gp["psi_F2"], psi_F2_true),
      delta_mu = diag$delta_mu, mean_bias = diag$mean_bias,
      r_kkT = diag$r_kkT, feas_gap = diag$feas_gap,
      fidelity = diag$fidelity, status = diag$status
    )
  }
  res_list <- list()

  df_true <- generate_sem_data(params$n, params$dist)
  miss_out <- apply_missingness_sem(df_true, params$miss, params$mech)
  df_miss <- miss_out$data
  act_m   <- miss_out$actual_miss

  # ── FIML Baseline ────────────────────────────────────────────────────────
  time_fiml <- system.time({
    fit_fiml <- tryCatch(lavaan::sem(sem_mod, data = df_miss, missing = "fiml"),
                         error = function(e) NULL)
    b_f <- NA; b_se_f <- NA; d_fiml <- NA
    gp <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    gs <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    if (!is.null(fit_fiml)) {
      bk <- extract_sem_path(fit_fiml)
      b_f <- bk["b"]; b_se_f <- bk["b_se"]
      implied_cov <- tryCatch(lavaan::lavInspect(fit_fiml, "cov.ov"),
                              error = function(e) NULL)
      if (!is.null(implied_cov)) d_fiml <- frob_dist(implied_cov, true_cov_sem)
      gp <- extract_sem_params(fit_fiml)
      gs <- extract_sem_params(fit_fiml, "se")
    }
  })["elapsed"]
  res_list[[1]] <- mk("FIML", d_fiml, b_f, b_se_f, time_fiml, gp, gs)

  # ── FIML + lavPredict: conditional-expectation completed data ─────────────
  time_fp <- system.time({
    b_fp <- NA; b_se_fp <- NA; d_fp <- NA
    gp <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    gs <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    fp_diag <- list()
    if (!is.null(fit_fiml)) {
      imp_fp <- tryCatch(lavaan::lavPredict(fit_fiml, type = "ov"),
                         error = function(e) NULL)
      if (!is.null(imp_fp)) {
        d_fp <- frob_dist(stats::cov(imp_fp), true_cov_sem)
        colnames(imp_fp) <- paste0("y", 1:6)
        fit_fp <- tryCatch(lavaan::sem(sem_mod, data = as.data.frame(imp_fp)),
                           error = function(e) NULL)
        bk <- extract_sem_path(fit_fp)
        b_fp <- bk["b"]; b_se_fp <- bk["b_se"]
        if (!is.null(fit_fp)) { gp <- extract_sem_params(fit_fp); gs <- extract_sem_params(fit_fp, "se") }
        fp_diag <- list(mean_bias = mean_bias_fn(imp_fp, df_true))
      }
    }
  })["elapsed"]
  res_list[[2]] <- mk("FIML_lavPredict", d_fp, b_fp, b_se_fp, time_fp,
                      gp, gs, diag = fp_diag)

  # ── MICE (random forest via ranger, MI m=20) ──────────────────────────────
  time_mice_rf <- system.time({
    m_mice_rf <- 20
    imp_mice_rf_list <- tryCatch({
      imp_obj_rf <- mice::mice(df_miss, m = m_mice_rf, method = "rf", printFlag = FALSE)
      mice::complete(imp_obj_rf, "all")
    }, error = function(e) NULL)

    b_mr <- NA; b_se_mr <- NA; d_mr <- NA
    gp <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    gps <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    mr_diag <- list()

    if (!is.null(imp_mice_rf_list)) {
      cov_list <- lapply(imp_mice_rf_list, function(x) stats::cov(x))
      avg_cov <- Reduce("+", cov_list) / length(cov_list)
      d_mr <- frob_dist(avg_cov, true_cov_sem)
      mr_diag <- list(mean_bias = mi_mean(mean_bias_fn, imp_mice_rf_list, df_true))

      fit_results <- lapply(imp_mice_rf_list, function(ds) {
        fit <- tryCatch(lavaan::sem(sem_mod, data = ds), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        list(est = extract_sem_params(fit, "est"), se = extract_sem_params(fit, "se"))
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
        b_mr <- gp["b21"]; b_se_mr <- gps["b21"]
      }
    }
  })["elapsed"]
  res_list[[3]] <- mk("MICE", d_mr, b_mr, b_se_mr, time_mice_rf,
                      gp, gps, diag = mr_diag)

  # ── missForest Baseline ───────────────────────────────────────────────────
  time_mf <- system.time({
    imp_mf <- tryCatch(missForest::missForest(df_miss, verbose = FALSE)$ximp,
                       error = function(e) NULL)
    b_mf <- NA; b_se_mf <- NA; d_mf <- NA
    gp <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    gs <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    mf_diag <- list()
    if (!is.null(imp_mf)) {
      d_mf <- frob_dist(stats::cov(imp_mf), true_cov_sem)
      fit_mf <- tryCatch(lavaan::sem(sem_mod, data = imp_mf), error = function(e) NULL)
      bk <- extract_sem_path(fit_mf)
      b_mf <- bk["b"]; b_se_mf <- bk["b_se"]
      if (!is.null(fit_mf)) { gp <- extract_sem_params(fit_mf); gs <- extract_sem_params(fit_mf, "se") }
      mf_diag <- list(mean_bias = mean_bias_fn(imp_mf, df_true))
    }
  })["elapsed"]
  res_list[[4]] <- mk("missForest", d_mf, b_mf, b_se_mf, time_mf,
                      gp, gs, diag = mf_diag)

  # ── prism_sem: FIML model-implied Σ target ───────────────────────────
  time_sf <- system.time({
    imp_sf <- tryCatch(suppressWarnings(
      # reuse the FIML fit from the baseline block when available to avoid
      # refitting the model per replication
      prism_sem(df_miss,
                model = if (!is.null(fit_fiml)) fit_fiml else sem_mod,
                initial_imputation = imp_mf, lambda_sigma = 1.0)),
      error = function(e) NULL)
    b_sf <- NA; b_se_sf <- NA; d_sf <- NA
    gp <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    gs <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    sf_diag <- list()
    if (!is.null(imp_sf)) {
      d_sf <- frob_dist(stats::cov(imp_sf[, paste0("y", 1:6)]), true_cov_sem)
      fit_sf <- tryCatch(lavaan::sem(sem_mod, data = imp_sf), error = function(e) NULL)
      bk <- extract_sem_path(fit_sf)
      b_sf <- bk["b"]; b_se_sf <- bk["b_se"]
      if (!is.null(fit_sf)) { gp <- extract_sem_params(fit_sf); gs <- extract_sem_params(fit_sf, "se") }
      sf_diag <- c(list(delta_mu = if (!is.null(imp_mf)) delta_mu_fn(imp_sf[, paste0("y", 1:6)], imp_mf) else NA_real_,
                        mean_bias = mean_bias_fn(imp_sf[, paste0("y", 1:6)], df_true)),
                   prism_diag_fn(imp_sf))
    }
  })["elapsed"]
  res_list[[5]] <- mk("PRISM", d_sf, b_sf, b_se_sf, time_sf, gp, gs,
                      pipeline_time = unname(time_mf) + unname(time_sf),
                      diag = sf_diag)

  # ── prism_mi: Proper Multiple Imputation ────────────────────────────────────
  time_smi <- system.time({
    m_prism <- 20
    imp_smi_list <- tryCatch(suppressWarnings({
      fit_fiml_base <- lavaan::sem(sem_mod, data = df_miss, missing = "fiml")
      if (two_level_mi) {
        prism_mi(df_miss, fit_fiml_base, m = m_prism, lambda_sigma = 1.0)
      } else {
        prism_mi(df_miss, fit_fiml_base, m = m_prism,
                 initial_imputation = imp_mf, lambda_sigma = 1.0)
      }
    }), error = function(e) NULL)

    b_smi <- NA; b_se_smi <- NA; d_smi <- NA
    gp <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    gps <- c(b21 = NA, psi_F1 = NA, psi_F2 = NA)
    smi_diag <- list()

    if (!is.null(imp_smi_list) && length(imp_smi_list) == m_prism) {
      cov_list <- lapply(imp_smi_list, function(x) stats::cov(x[, paste0("y", 1:6)]))
      avg_cov <- Reduce("+", cov_list) / length(cov_list)
      d_smi <- frob_dist(avg_cov, true_cov_sem)
      smi_diag <- c(list(delta_mu  = mi_mean(delta_mu_fn,  imp_smi_list, imp_mf),
                         mean_bias = mi_mean(mean_bias_fn, imp_smi_list, df_true)),
                    mi_prism_diag_fn(imp_smi_list))

      fit_results <- lapply(imp_smi_list, function(ds) {
        fit <- tryCatch(lavaan::sem(sem_mod, data = ds), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        list(est = extract_sem_params(fit, "est"), se = extract_sem_params(fit, "se"))
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
        b_smi <- gp["b21"]; b_se_smi <- gps["b21"]
      }
    }
  })["elapsed"]
  res_list[[6]] <- mk("PRISM_MI", d_smi, b_smi, b_se_smi, time_smi,
                      gp, gps,
                      pipeline_time = unname(time_mf) + unname(time_smi),
                      diag = smi_diag)

  rm(df_true, df_miss, imp_mice_rf_list, imp_mf, imp_sf, imp_fp, imp_smi_list)
  do.call(rbind, res_list)
}

# ── Worker Dispatcher ─────────────────────────────────────────────────────────
run_iteration <- function(sim_id, params) {
  if (identical(params$model, "SEM")) return(run_iteration_sem(sim_id, params))
  run_iteration_gcm(sim_id, params)
}

# ── SLURM Array Dispatch ─────────────────────────────────────────────────────
conditions <- expand.grid(model = grid_model, n = grid_n, miss = grid_miss,
                          dist = grid_dist, mech = grid_mech,
                          stringsAsFactors = FALSE)
total_conditions <- nrow(conditions)

array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.na(array_id) && array_id >= 1 && array_id <= total_conditions) {
  current_conditions <- conditions[array_id, , drop = FALSE]
  output_file <- sprintf("sim_results/prod_results_%d.rds", array_id)
} else if (!is.na(array_id)) {
  stop(sprintf(
    "SLURM_ARRAY_TASK_ID = %d is outside the grid of %d conditions. ",
    "Update the --array directive in simulations/submit_simulation.slurm ",
    "to match the condition grid in simulations/simulation.R.",
    array_id, total_conditions
  ))
} else {
  current_conditions <- conditions
  output_file <- "sim_results/prod_results.rds"
}

# ── Execution ────────────────────────────────────────────────────────────────
cat(sprintf("Grid: %d conditions (%d models) × %d reps × 6 methods = %d total rows\n",
            total_conditions, length(grid_model), n_sims,
            total_conditions * n_sims * 6))
cat(sprintf("Parallel cores: %d\n", num_cores))
cat(sprintf("Output: %s\n\n", output_file))

dir.create("sim_results", showWarnings = FALSE)
dir.create("sim_raw_data", showWarnings = FALSE)

all_results <- list()

for (i in seq_len(nrow(current_conditions))) {
  params <- current_conditions[i, ]
  cond_idx <- if (!is.na(array_id)) array_id else i
  cat(sprintf("[%s] Condition %d/%d: Model=%s, N=%d, Miss=%.2f, Dist=%s, Mech=%s\n",
              Sys.time(), cond_idx, total_conditions, params$model, params$n,
              params$miss, params$dist, params$mech))

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
