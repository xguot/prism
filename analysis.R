# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Analysis for prism (v2, GCM + SEM)
# Style adapted from mda/analysis.R — parameter-level breakdown, MSE, heatmaps
# Reads the per-task files sim_results/prod_results_<id>.rds via load_results.R
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)

source("load_results.R")
prod <- load_results()

# ── Shared helpers ────────────────────────────────────────────────────────────
hr <- function(label) {
  cat("\n", strrep("=", 90), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("=", 90), "\n\n", sep = "")
}

# Pooled summaries use the median and MAD (SD-consistent) instead of the
# mean and SD: in the heavy-tail cells (SEM x Outlier / Lognormal) a small
# number of degenerate fits produces estimates many orders of magnitude
# from the truth, and mean-based pooling is dominated by them.  Median/MAD
# describe the bulk of the replications.  Counts and timing keep their
# original (non-robust) definitions.
compute_metrics <- function(results, param, true_val,
                           ec_map = est_cols, sc_map = se_cols) {
  ec <- ec_map[param]
  sc <- sc_map[param]
  has_se <- sc %in% names(results)

  valid <- results %>% filter(converged == 1, is.finite(.data[[ec]]))
  if (has_se) valid <- valid %>% filter(is.finite(.data[[sc]]))

  out <- valid %>%
    group_by(method) %>%
    summarise(
      n_converged      = n(),
      median_estimate  = stats::median(.data[[ec]], na.rm = TRUE),
      mad_est          = stats::mad(.data[[ec]], na.rm = TRUE),
      avg_model_se     = if (has_se) mean(.data[[sc]], na.rm = TRUE) else NA_real_,
      .groups          = "drop"
    ) %>%
    mutate(
      bias     = median_estimate - true_val,
      se_ratio = if (has_se) avg_model_se / mad_est else NA_real_,
      param    = param
    )

  out
}

# GCM truth: sigma2_i = 1, sigma2_s = 1, sigma_is = 0, mu_i = 6, mu_s = 2
beta_true <- c(
  var_intercept = 1, var_slope = 1, cov_intercept_slope = 0,
  mean_intercept = 6, mean_slope = 2
)
param_names  <- names(beta_true)
param_labels <- c(
  var_intercept       = "Variance (Intercept)",
  var_slope           = "Variance (Slope)",
  cov_intercept_slope = "Covariance (Int, Slope)",
  mean_intercept      = "Mean (Intercept)",
  mean_slope          = "Mean (Slope)"
)

# SEM truth: b21 = 0.5, psi_F1 = 1, psi_F2 = 0.75 (residual variance of eta2)
sem_true <- c(b21 = 0.5, psi_F1 = 1, psi_F2 = 0.75)
sem_names  <- names(sem_true)
sem_labels <- c(
  b21    = "Path eta2 ~ eta1",
  psi_F1 = "Variance (Factor 1)",
  psi_F2 = "Residual Variance (Factor 2)"
)

# Map production column names to GCM parameter names
est_cols <- c(
  var_intercept  = "est_var_L",  var_slope = "est_var_S",
  cov_intercept_slope = "est_cov_LS",
  mean_intercept = "est_L",      mean_slope = "est_S"
)
se_cols <- c(
  var_intercept  = "se_var_L",  var_slope = "se_var_S",
  cov_intercept_slope = "se_cov_LS",
  mean_intercept = "se_L",      mean_slope = "se_S"
)

# Map production column names to SEM parameter names
sem_est_cols <- c(b21 = "est_b21", psi_F1 = "est_var_F1", psi_F2 = "est_var_F2")
sem_se_cols  <- c(b21 = "se_b21",  psi_F1 = "se_var_F1",  psi_F2 = "se_var_F2")

cat(sprintf("Production data: %d rows, %d methods, %d models\n",
            nrow(prod), length(unique(prod$method)), length(unique(prod$model))))

# ── Build per-parameter long tables (per model) ───────────────────────────────
gcm_prod <- prod %>% filter(model == "GCM")
sem_prod <- prod %>% filter(model == "SEM")

all_params <- do.call(rbind, lapply(param_names, function(pn) {
  tv <- beta_true[pn]
  ec <- est_cols[pn]
  gcm_prod %>%
    mutate(
      param    = pn,
      est      = .data[[ec]],
      bias_raw = est - tv,
      relbias  = if (abs(tv) < 1e-12) est - tv else 100 * (est - tv) / tv
    ) %>%
    select(sim_id, N, miss, actual_miss, dist, mech, method, param, est, bias_raw, relbias)
}))

sem_all_params <- do.call(rbind, lapply(sem_names, function(pn) {
  tv <- sem_true[pn]
  ec <- sem_est_cols[pn]
  sem_prod %>%
    mutate(
      param    = pn,
      est      = .data[[ec]],
      bias_raw = est - tv,
      relbias  = 100 * (est - tv) / tv
    ) %>%
    select(sim_id, N, miss, actual_miss, dist, mech, method, param, est, bias_raw, relbias)
}))

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 0 — Convergence, Empirical SE, Model SE, SE Ratio (all parameters)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 0G — GCM: Median Estimate, MAD, Model SE, SE Ratio")

for (pn in param_names) {
  cat(sprintf("\n--- %s (truth = %.0f) ---\n", pn, beta_true[pn]))
  sm <- compute_metrics(gcm_prod, pn, beta_true[pn])
  print(as.data.frame(sm), row.names = FALSE)
}

hr("TABLE 0S — SEM: Median Estimate, MAD, Model SE, SE Ratio")

for (pn in sem_names) {
  cat(sprintf("\n--- %s (truth = %.2f) ---\n", pn, sem_true[pn]))
  sm <- compute_metrics(sem_prod, pn, sem_true[pn],
                        ec_map = sem_est_cols, sc_map = sem_se_cols)
  print(as.data.frame(sm), row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 1 — Frobenius Distance (primary metric, per model)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 1 — Frobenius Distance to True Covariance (lower = better)")

prod %>%
  group_by(model, dist, mech, method) %>%
  summarise(
    Frobenius = stats::median(f_dist, na.rm = TRUE),
    MAD       = stats::mad(f_dist, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f (%.2f)", Frobenius, MAD)) %>%
  select(model, dist, mech, method, Display) %>%
  pivot_wider(names_from = mech, values_from = Display) %>%
  arrange(model, dist, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2 — Frobenius by Sample Size (MAR only)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 2 — Frobenius Distance by Sample Size (MAR only)")

prod %>%
  filter(mech == "MAR") %>%
  group_by(model, N, dist, method) %>%
  summarise(Frobenius = stats::median(f_dist, na.rm = TRUE), .groups = "drop") %>%
  mutate(Display = sprintf("%.2f", Frobenius)) %>%
  select(model, N, dist, method, Display) %>%
  pivot_wider(names_from = N, values_from = Display, names_prefix = "N=") %>%
  arrange(model, dist, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 3 — Frobenius by Missingness Rate (MAR only)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 3 — Frobenius Distance by Missingness Rate (MAR only)")

prod %>%
  filter(mech == "MAR") %>%
  mutate(miss_pct = sprintf("%.0f%%", miss * 100)) %>%
  group_by(model, miss_pct, dist, method) %>%
  summarise(Frobenius = stats::median(f_dist, na.rm = TRUE), .groups = "drop") %>%
  mutate(Display = sprintf("%.2f", Frobenius)) %>%
  select(model, miss_pct, dist, method, Display) %>%
  pivot_wider(names_from = miss_pct, values_from = Display) %>%
  arrange(model, dist, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 4 — Relative Bias per Parameter (MAR, aggregated across N x miss)
# ══════════════════════════════════════════════════════════════════════════════
for (pn in param_names) {
  tv <- beta_true[pn]
  title <- if (abs(tv) < 1e-12)
    sprintf("TABLE 4%s — Raw Bias (median (MAD)): %s (truth = 0)", pn, param_labels[pn])
  else
    sprintf("TABLE 4%s — Relative Bias (median (MAD)): %s (%%, truth = %.0f)", pn, param_labels[pn], tv)
  hr(title)

  all_params %>%
    filter(param == pn, mech == "MAR") %>%
    group_by(dist, method) %>%
    summarise(
      M   = stats::median(relbias, na.rm = TRUE),
      MAD = stats::mad(relbias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Display = sprintf("%+.2f (%.2f)", M, MAD)) %>%
    select(dist, method, Display) %>%
    pivot_wider(names_from = dist, values_from = Display) %>%
    arrange(method) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

for (pn in sem_names) {
  tv <- sem_true[pn]
  hr(sprintf("TABLE 4%s — Relative Bias (median (MAD)): %s (%%, truth = %.2f)", pn, sem_labels[pn], tv))

  sem_all_params %>%
    filter(param == pn, mech == "MAR") %>%
    group_by(dist, method) %>%
    summarise(
      M   = stats::median(relbias, na.rm = TRUE),
      MAD = stats::mad(relbias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Display = sprintf("%+.2f (%.2f)", M, MAD)) %>%
    select(dist, method, Display) %>%
    pivot_wider(names_from = dist, values_from = Display) %>%
    arrange(method) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 5 — MSE per Parameter (MAR, aggregated)
# ══════════════════════════════════════════════════════════════════════════════
for (pn in param_names) {
  tv <- beta_true[pn]
  hr(sprintf("TABLE 5%s — Median-based MSE: %s (truth = %.0f)", pn, param_labels[pn], tv))

  all_params %>%
    filter(param == pn, mech == "MAR") %>%
    group_by(dist, method) %>%
    summarise(
      bias_med = stats::median(est, na.rm = TRUE) - tv,
      MAD      = stats::mad(est, na.rm = TRUE),
      MSE      = bias_med^2 + MAD^2,
      .groups  = "drop"
    ) %>%
    mutate(Display = sprintf("%.4f", MSE)) %>%
    select(dist, method, Display) %>%
    pivot_wider(names_from = dist, values_from = Display) %>%
    arrange(method) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

for (pn in sem_names) {
  tv <- sem_true[pn]
  hr(sprintf("TABLE 5%s — Median-based MSE: %s (truth = %.2f)", pn, sem_labels[pn], tv))

  sem_all_params %>%
    filter(param == pn, mech == "MAR") %>%
    group_by(dist, method) %>%
    summarise(
      bias_med = stats::median(est, na.rm = TRUE) - tv,
      MAD      = stats::mad(est, na.rm = TRUE),
      MSE      = bias_med^2 + MAD^2,
      .groups  = "drop"
    ) %>%
    mutate(Display = sprintf("%.4f", MSE)) %>%
    select(dist, method, Display) %>%
    pivot_wider(names_from = dist, values_from = Display) %>%
    arrange(method) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 6 — Outlier Impact: Frobenius Degradation (MAR)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 6 — Outlier Impact: Delta Frobenius (Normal -> Outlier, MAR)")

prod %>%
  filter(mech == "MAR", dist %in% c("Normal", "Outlier")) %>%
  group_by(model, dist, method) %>%
  summarise(Frob = stats::median(f_dist, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = dist, values_from = Frob) %>%
  mutate(
    Delta   = Outlier - Normal,
    Display = sprintf("%.2f -> %.2f  (Delta %+.2f)", Normal, Outlier, Delta)
  ) %>%
  arrange(model, Delta) %>%
  select(model, method, Display) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 7 — Timing Summary (MAR, mean ± SD seconds)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 7 — Timing Summary (mean +/- SD seconds)")

prod %>%
  filter(mech == "MAR") %>%
  group_by(model, method) %>%
  summarise(
    Mean = mean(time_sec, na.rm = TRUE),
    SD   = sd(time_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f +/- %.2f", Mean, SD)) %>%
  select(model, method, Display) %>%
  arrange(model, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 8 — Relative Bias Heatmap Summary (MAR, pooled across conditions)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 8G — Median Relative Bias Summary (GCM, MAR, pooled across N x miss)")

all_params %>%
  filter(mech == "MAR") %>%
  group_by(param, method) %>%
  summarise(RelBias = stats::median(relbias, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = RelBias) %>%
  mutate(across(where(is.numeric), ~ sprintf("%+.1f%%", .x))) %>%
  arrange(method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

hr("TABLE 8S — Median Relative Bias Summary (SEM, MAR, pooled across N x miss)")

sem_all_params %>%
  filter(mech == "MAR") %>%
  group_by(param, method) %>%
  summarise(RelBias = stats::median(relbias, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = RelBias) %>%
  mutate(across(where(is.numeric), ~ sprintf("%+.1f%%", .x))) %>%
  arrange(method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 9 — PRISM Diagnostics and Mean Drift (MAR, pooled across conditions)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 9 — PRISM Diagnostics: KKT Residual, Feasibility Gap, Fidelity, Certification (MAR)")

prod %>%
  filter(mech == "MAR", method %in% c("PRISM", "PRISM_MI")) %>%
  group_by(model, method, dist) %>%
  summarise(
    r_kkT    = stats::median(r_kkT, na.rm = TRUE),
    feas_gap = stats::median(feas_gap, na.rm = TRUE),
    fidelity = stats::median(fidelity, na.rm = TRUE),
    frac_converged = mean(status == "converged", na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    r_kkT          = sprintf("%.2e", r_kkT),
    feas_gap       = sprintf("%.3f", feas_gap),
    fidelity       = sprintf("%.3f", fidelity),
    frac_converged = sprintf("%.2f", frac_converged)
  ) %>%
  arrange(model, method, dist) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 10 — Mean Drift and Mean Bias (MAR, pooled across conditions)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 10 — Mean Drift (delta_mu) and Mean Bias (MAR)")

prod %>%
  filter(mech == "MAR") %>%
  group_by(model, method) %>%
  summarise(
    delta_mu  = stats::median(delta_mu, na.rm = TRUE),
    mean_bias = stats::median(mean_bias, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    delta_mu  = sprintf("%.2e", delta_mu),
    mean_bias = sprintf("%.4f", mean_bias)
  ) %>%
  arrange(model, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 11 — PRISM Stationarity Certification (MAR)
# ══════════════════════════════════════════════════════════════════════════════
# Fraction of PRISM replications that reached first-order stationarity
# (status = "converged") rather than stopping at max_iter.  The headline
# tables above are NOT filtered on status; this table shows how much of the
# evidence is certified in each cell.
hr("TABLE 11 — PRISM Stationarity Certification (MAR, frac converged)")

prod %>%
  filter(mech == "MAR", method %in% c("PRISM", "PRISM_MI")) %>%
  mutate(miss_pct = sprintf("%.0f%%", miss * 100)) %>%
  group_by(model, method, dist, miss_pct) %>%
  summarise(frac = mean(status == "converged", na.rm = TRUE), .groups = "drop") %>%
  mutate(Display = sprintf("%.2f", frac)) %>%
  select(model, method, dist, miss_pct, Display) %>%
  pivot_wider(names_from = miss_pct, values_from = Display) %>%
  arrange(model, method, dist) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
