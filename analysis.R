# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Analysis for prism
# Style adapted from mda/analysis.R — parameter-level breakdown, MSE, heatmaps
# Uses prod_results.rds (post-HPC)
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)

prod  <- readRDS("sim_results/prod_results.rds")

# ── Shared helpers ────────────────────────────────────────────────────────────
hr <- function(label) {
  cat("\n", strrep("=", 90), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("=", 90), "\n\n", sep = "")
}

compute_metrics <- function(results, param, true_val) {
  ec <- est_cols[param]
  sc <- se_cols[param]
  has_se <- sc %in% names(results)

  valid <- results %>% filter(converged == 1, is.finite(.data[[ec]]))
  if (has_se) valid <- valid %>% filter(is.finite(.data[[sc]]))

  out <- valid %>%
    group_by(method) %>%
    summarise(
      n_converged   = n(),
      mean_estimate = mean(.data[[ec]], na.rm = TRUE),
      empirical_se  = sd(.data[[ec]], na.rm = TRUE),
      avg_model_se  = if (has_se) mean(.data[[sc]], na.rm = TRUE) else NA_real_,
      .groups       = "drop"
    ) %>%
    mutate(
      bias     = mean_estimate - true_val,
      se_ratio = if (has_se) avg_model_se / empirical_se else NA_real_,
      param    = param
    )

  out
}
beta_true   <- c(var_intercept = 1, var_slope = 1, cov_intercept_slope = 0, mean_intercept = 6, mean_slope = 2)
param_names <- names(beta_true)
param_labels <- c(
  var_intercept  = "Variance (Intercept)",
  var_slope  = "Variance (Slope)",
  cov_intercept_slope = "Covariance (Int, Slope)",
  mean_intercept = "Mean (Intercept)",
  mean_slope = "Mean (Slope)"
)

cat(sprintf("Production data: %d rows, %d methods\n",
            nrow(prod), length(unique(prod$method))))

# ── Build per-parameter long table (all methods, all conditions) ─────────────
# Map production column names to GCM parameter names
est_cols <- c(
  var_intercept  = "est_var_L",  var_slope  = "est_var_S", cov_intercept_slope = "est_cov_LS",
  mean_intercept = "est_L",       mean_slope = "est_S"
)
se_cols <- c(
  var_intercept  = "se_var_L",  var_slope  = "se_var_S", cov_intercept_slope = "se_cov_LS",
  mean_intercept = "se_L",       mean_slope = "se_S"
)

all_params <- do.call(rbind, lapply(param_names, function(pn) {
  tv <- beta_true[pn]
  ec <- est_cols[pn]
  prod %>%
    mutate(
      param    = pn,
      est      = .data[[ec]],
      bias_raw = est - tv,
      relbias  = if (abs(tv) < 1e-12) est - tv else 100 * (est - tv) / tv
    ) %>%
    select(sim_id, N, miss, dist, mech, method, param, est, bias_raw, relbias)
}))

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 0 — Convergence, Empirical SE, Model SE, SE Ratio (all parameters)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 0 — Convergence, Empirical SE, Model SE, SE Ratio (all parameters)")

for (pn in param_names) {
  cat(sprintf("\n--- %s (truth = %.0f) ---\n", pn, beta_true[pn]))
  sm <- compute_metrics(prod, pn, beta_true[pn])
  print(as.data.frame(sm), row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 1 — Frobenius Distance (primary metric)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 1 — Frobenius Distance to True Covariance (lower = better)")

prod %>%
  group_by(dist, mech, method) %>%
  summarise(
    Frobenius = mean(f_dist, na.rm = TRUE),
    SD        = sd(f_dist, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f (%.2f)", Frobenius, SD)) %>%
  select(dist, mech, method, Display) %>%
  pivot_wider(names_from = mech, values_from = Display) %>%
  arrange(dist, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2 — Frobenius by Sample Size (MAR only)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 2 — Frobenius Distance by Sample Size (MAR only)")

prod %>%
  filter(mech == "MAR") %>%
  group_by(N, dist, method) %>%
  summarise(Frobenius = mean(f_dist, na.rm = TRUE), .groups = "drop") %>%
  mutate(Display = sprintf("%.2f", Frobenius)) %>%
  select(N, dist, method, Display) %>%
  pivot_wider(names_from = N, values_from = Display, names_prefix = "N=") %>%
  arrange(dist, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 3 — Frobenius by Missingness Rate (MAR only)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 3 — Frobenius Distance by Missingness Rate (MAR only)")

prod %>%
  filter(mech == "MAR") %>%
  mutate(miss_pct = sprintf("%.0f%%", miss * 100)) %>%
  group_by(miss_pct, dist, method) %>%
  summarise(Frobenius = mean(f_dist, na.rm = TRUE), .groups = "drop") %>%
  mutate(Display = sprintf("%.2f", Frobenius)) %>%
  select(miss_pct, dist, method, Display) %>%
  pivot_wider(names_from = miss_pct, values_from = Display) %>%
  arrange(dist, method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 4 — Relative Bias per GCM Parameter (MAR, aggregated across N × miss)
# ══════════════════════════════════════════════════════════════════════════════
for (pn in param_names) {
  tv <- beta_true[pn]
  title <- if (abs(tv) < 1e-12)
    sprintf("TABLE 4%s — Raw Bias: %s (truth = 0)", pn, param_labels[pn])
  else
    sprintf("TABLE 4%s — Relative Bias: %s (%%, truth = %.0f)", pn, param_labels[pn], tv)
  hr(title)

  all_params %>%
    filter(param == pn, mech == "MAR") %>%
    group_by(dist, method) %>%
    summarise(
      M  = mean(relbias, na.rm = TRUE),
      SD = sd(relbias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Display = sprintf("%+.2f (%.2f)", M, SD)) %>%
    select(dist, method, Display) %>%
    pivot_wider(names_from = dist, values_from = Display) %>%
    arrange(method) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 5 — MSE per GCM Parameter (MAR, aggregated)
# ══════════════════════════════════════════════════════════════════════════════
for (pn in param_names) {
  tv <- beta_true[pn]
  hr(sprintf("TABLE 5%s — MSE: %s (truth = %.0f)", pn, param_labels[pn], tv))

  all_params %>%
    filter(param == pn, mech == "MAR") %>%
    group_by(dist, method) %>%
    summarise(
      bias_raw = mean(est, na.rm = TRUE) - tv,
      ESE      = sd(est, na.rm = TRUE),
      MSE      = bias_raw^2 + ESE^2,
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
  group_by(dist, method) %>%
  summarise(Frob = mean(f_dist, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = dist, values_from = Frob) %>%
  mutate(
    Delta   = Outlier - Normal,
    Display = sprintf("%.2f -> %.2f  (Delta %+.2f)", Normal, Outlier, Delta)
  ) %>%
  arrange(Delta) %>%
  select(method, Display) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 7 — Timing Summary (MAR, mean ± SD seconds)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 7 — Timing Summary (mean +/- SD seconds)")

prod %>%
  filter(mech == "MAR") %>%
  group_by(method) %>%
  summarise(
    Mean = mean(time_sec, na.rm = TRUE),
    SD   = sd(time_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f +/- %.2f", Mean, SD)) %>%
  select(method, Display) %>%
  arrange(method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 8 — Relative Bias Heatmap Summary (MAR, pooled across conditions)
# ══════════════════════════════════════════════════════════════════════════════
hr("TABLE 8 — Relative Bias Summary (MAR, pooled across N x miss)")

all_params %>%
  filter(mech == "MAR") %>%
  group_by(param, method) %>%
  summarise(RelBias = mean(relbias, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = param, values_from = RelBias) %>%
  mutate(across(where(is.numeric), ~ sprintf("%+.1f%%", .x))) %>%
  arrange(method) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
