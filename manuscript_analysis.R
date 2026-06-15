# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Analysis for smriti
# Style adapted from mda/analysis.R — parameter-level breakdown, MSE, heatmaps
# Uses prod_results.rds (post-HPC) and tune_results.rds
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)

prod  <- readRDS("sim_results/prod_results.rds")
tune_file <- "sim_results/tune_results.rds"
tune_exists <- file.exists(tune_file)

# ── Shared helpers ────────────────────────────────────────────────────────────
hr <- function(label) {
  cat("\n", strrep("=", 90), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("=", 90), "\n\n", sep = "")
}

beta_true   <- c(psi_L = 1, psi_S = 1, psi_LS = 0, beta_L = 6, beta_S = 2)
param_names <- names(beta_true)
param_labels <- c(
  psi_L  = expression(psi[L] ~ "(Var intercept)"),
  psi_S  = expression(psi[S] ~ "(Var slope)"),
  psi_LS = expression(psi[LS] ~ "(Covariance)"),
  beta_L = expression(beta[L] ~ "(Intercept)"),
  beta_S = expression(beta[S] ~ "(Slope)")
)

cat(sprintf("Production data: %d rows, %d methods\n",
            nrow(prod), length(unique(prod$method))))
if (tune_exists) {
  tune <- readRDS(tune_file)
  cat(sprintf("Tuning data:     %d rows, lambda in {%s}\n",
              nrow(tune), paste(sort(unique(tune$lambda)), collapse = ", ")))
} else {
  tune <- NULL
  warning("tune_results.rds not found. Skipping tuning-dependent tables.")
}

# ── Build per-parameter long table (all methods, all conditions) ─────────────
# Map production column names to GCM parameter names
est_cols <- c(
  psi_L  = "est_var_L",  psi_S  = "est_var_S", psi_LS = "est_cov_LS",
  beta_L = "est_L",       beta_S = "est_S"
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

# ══════════════════════════════════════════════════════════════════════════════
# Tuning-dependent tables
# ══════════════════════════════════════════════════════════════════════════════
if (tune_exists) {
  hr("TABLE T1 — Lambda Tuning: Frobenius by lambda x robust (MAR, pooled)")

  tune %>%
    filter(mech == "MAR") %>%
    group_by(lambda, robust, dist) %>%
    summarise(
      Frobenius = mean(f_dist, na.rm = TRUE),
      Bias_pct  = mean(rel_bias, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      Robust  = ifelse(robust, "Robust", "Pearson"),
      Display = sprintf("lambda=%.2f %-7s  Frob=%.2f  Bias=%+.1f%%",
                        lambda, Robust, Frobenius, Bias_pct)
    ) %>%
    arrange(dist, lambda, desc(robust)) -> tune_tbl

  for (d in unique(tune_tbl$dist)) {
    cat(sprintf("\n  -- %s --\n", d))
    subset_rows <- tune_tbl %>% filter(dist == d) %>% pull(Display)
    for (line in subset_rows) cat("  ", line, "\n")
  }

  hr("TABLE T2 — FIML vs Smriti_FIML (Tuning Study Summary)")
  tune %>%
    filter(mech == "MAR", method %in% c("FIML", "Smriti_FIML")) %>%
    group_by(method) %>%
    summarise(
      Frob = mean(f_dist, na.rm = TRUE),
      SD   = sd(f_dist, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}
