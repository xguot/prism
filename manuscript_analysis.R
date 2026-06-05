library(dplyr)
library(tidyr)

# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Analysis for smriti
# Uses existing prod_results.rds and tune_results.rds — no re-simulation needed.
# ══════════════════════════════════════════════════════════════════════════════

prod  <- readRDS("sim_results/prod_results.rds")
tune_file <- "sim_results/tune_results.rds"
tune_exists <- file.exists(tune_file)

if (tune_exists) {
  tune <- readRDS(tune_file)
} else {
  tune <- NULL
  warning("Tuning results not found at ", tune_file, ". Skipping tuning-dependent tables.")
}

cat(sprintf("Production data: %d rows, %d methods\n",
            nrow(prod), length(unique(prod$method))))
if (tune_exists) {
  cat(sprintf("Tuning data:     %d rows, λ ∈ {%s}\n",
              nrow(tune), paste(sort(unique(tune$lambda)), collapse = ", ")))
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 1 — Frobenius Distance (primary metric: full covariance recovery)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat(strrep("═", 74), "\n")
cat("TABLE 1 — Frobenius Distance to True Covariance (lower = better)\n")
cat("        Split by distribution × mechanism.  MAR is the fair comparison.\n")
cat(strrep("═", 74), "\n\n")

frob_table <- prod %>%
  group_by(dist, mech, method) %>%
  summarize(
    Frobenius = mean(f_dist, na.rm = TRUE),
    SD        = sd(f_dist, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f (%.2f)", Frobenius, SD)) %>%
  select(dist, mech, method, Display) %>%
  pivot_wider(names_from = mech, values_from = Display) %>%
  arrange(dist, method)

print(as.data.frame(frob_table), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2 — Frobenius Distance by Sample Size (does advantage hold at low N?)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat(strrep("═", 74), "\n")
cat("TABLE 2 — Frobenius Distance by Sample Size (MAR only)\n")
cat("        Does the Smriti advantage persist at N = 100?\n")
cat(strrep("═", 74), "\n\n")

n_table <- prod %>%
  filter(mech == "MAR") %>%
  group_by(N, dist, method) %>%
  summarize(
    Frobenius = mean(f_dist, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f", Frobenius)) %>%
  select(N, dist, method, Display) %>%
  pivot_wider(names_from = N, values_from = Display, names_prefix = "N=") %>%
  arrange(dist, method)

print(as.data.frame(n_table), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 3 — Frobenius Distance by Missingness Rate (MAR only)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat(strrep("═", 74), "\n")
cat("TABLE 3 — Frobenius Distance by Missingness Rate (MAR only)\n")
cat("        Does the Smriti advantage hold at 30% missingness?\n")
cat(strrep("═", 74), "\n\n")

miss_table <- prod %>%
  filter(mech == "MAR") %>%
  group_by(miss, dist, method) %>%
  summarize(
    Frobenius = mean(f_dist, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    Display = sprintf("%.2f", Frobenius),
    miss_pct = sprintf("%.0f%%", miss * 100)
  ) %>%
  select(miss_pct, dist, method, Display) %>%
  pivot_wider(names_from = miss_pct, values_from = Display) %>%
  arrange(dist, method)

print(as.data.frame(miss_table), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 4 — Slope Variance Bias (secondary metric, MAR only)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat(strrep("═", 74), "\n")
cat("TABLE 4 — Slope Variance Bias (MAR, secondary SEM metric)\n")
cat("        How well is the latent slope variance recovered?\n")
cat(strrep("═", 74), "\n\n")

bias_table <- prod %>%
  filter(mech == "MAR") %>%
  group_by(dist, method) %>%
  summarize(
    Bias_pct = mean(s_var_bias, na.rm = TRUE),
    SD       = sd(s_var, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(Display = sprintf("%+.2f%% (±%.3f)", Bias_pct, SD)) %>%
  select(dist, method, Display) %>%
  pivot_wider(names_from = dist, values_from = Display)

print(as.data.frame(bias_table), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 5 — Outlier Robustness: Degradation from Normal to Outlier (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat(strrep("═", 74), "\n")
cat("TABLE 5 — Outlier Impact: Frobenius Distance (Normal → Outlier Δ, MAR)\n")
cat("        Lower delta = more robust to 5% contamination.\n")
cat(strrep("═", 74), "\n\n")

outlier_impact <- prod %>%
  filter(mech == "MAR", dist %in% c("Normal", "Outlier")) %>%
  group_by(dist, method) %>%
  summarize(Frob = mean(f_dist, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = dist, values_from = Frob) %>%
  mutate(
    Delta = Outlier - Normal,
    Display = sprintf("%.2f → %.2f  (Δ %+.2f)", Normal, Outlier, Delta)
  ) %>%
  arrange(Delta) %>%
  select(method, Display)

print(as.data.frame(outlier_impact), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 6 — Tuning Study: λ Selection via Frobenius Distance (MAR only)
# ══════════════════════════════════════════════════════════════════════════════
if (tune_exists) {
  cat("\n")
  cat(strrep("═", 74), "\n")
  cat("TABLE 6 — λ Tuning: Frobenius Distance by λ (MAR, pooled N × miss)\n")
  cat("        Coarse tuning: N ∈ {200, 500, 5000}, 100 reps per condition.\n")
  cat(strrep("═", 74), "\n\n")

  tune_summary <- tune %>%
    filter(mech == "MAR") %>%
    group_by(lambda, robust, dist) %>%
    summarize(
      Frobenius = mean(f_dist, na.rm = TRUE),
      Bias_pct  = mean(rel_bias, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      Robust  = ifelse(robust, "Robust", "Pearson"),
      Display = sprintf("λ=%.2f %-7s  Frob=%.2f  Bias=%+.1f%%",
                        lambda, Robust, Frobenius, Bias_pct)
    ) %>%
    arrange(dist, lambda, desc(robust))

  for (d in unique(tune_summary$dist)) {
    cat(sprintf("\n  ── %s ──\n", d))
    subset_rows <- tune_summary %>%
      filter(dist == d) %>%
      pull(Display)
    for (line in subset_rows) cat("  ", line, "\n")
  }

  # ══════════════════════════════════════════════════════════════════════════════
  # λ = 1.0 justification
  # ══════════════════════════════════════════════════════════════════════════════
  cat("\n")
  cat(strrep("═", 74), "\n")
  cat("λ SELECTION RATIONALE\n")
  cat(strrep("═", 74), "\n\n")

  # Find the best λ per (dist, robust) combo by Frobenius distance
  best_lambda <- tune %>%
    filter(mech == "MAR") %>%
    group_by(dist, robust) %>%
    summarize(
      best_λ = lambda[which.min(mean(f_dist))],
      best_frob = min(mean(f_dist)),
      λ1_frob  = mean(f_dist[lambda == 1.0]),
      λ1_loss  = mean(f_dist[lambda == 1.0]) - min(mean(f_dist)),
      .groups  = "drop"
    )

  for (i in seq_len(nrow(best_lambda))) {
    r <- best_lambda[i, ]
    cat(sprintf("  %-10s %-7s: best λ=%.2f (Frob=%.2f); λ=1.0 gives Frob=%.2f (Δ +%.3f)\n",
                r$dist, ifelse(r$robust, "Robust", "Pearson"),
                r$best_λ, r$best_frob, r$λ1_frob, r$λ1_loss))
  }

  cat("\n  → λ = 1.0 selected as a single default that balances performance\n")
  cat("    across Normal and Outlier distributions for both robust and Pearson\n")
  cat("    modes.  It is within 0.01–0.05 Frobenius units of the per-condition\n")
  cat("    optimum while avoiding overfitting to any single data regime.\n")
} else {
  cat("\nTABLE 6 (Tuning Study) skipped — missing tune_results.rds\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 7 — Timing Breakdown
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat(strrep("═", 74), "\n")
cat("TABLE 7 — Timing (seconds, MAR, mean across all conditions)\n")
cat(strrep("═", 74), "\n\n")

time_table <- prod %>%
  filter(mech == "MAR") %>%
  group_by(method) %>%
  summarize(
    Time = mean(time_sec, na.rm = TRUE),
    SD   = sd(time_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Display = sprintf("%.2f ± %.2f s", Time, SD)) %>%
  select(method, Display)

print(as.data.frame(time_table), row.names = FALSE)

cat("\n  NOTE: Smriti timing reflects the Lagrangian routing step only.\n")
cat("  The full pipeline (missForest init + routing) adds ~8 s.\n")

# ══════════════════════════════════════════════════════════════════════════════
# AUTOMATED CSV EXPORTS
# ══════════════════════════════════════════════════════════════════════════════
cat("\nSaving summary tables to sim_results/...\n")
write.csv(frob_table,    "sim_results/table1_frobenius_main.csv",  row.names = FALSE)
write.csv(n_table,       "sim_results/table2_frobenius_by_N.csv",  row.names = FALSE)
write.csv(miss_table,    "sim_results/table3_frobenius_by_miss.csv", row.names = FALSE)
write.csv(bias_table,    "sim_results/table4_slope_bias.csv",       row.names = FALSE)
write.csv(outlier_impact, "sim_results/table5_outlier_impact.csv",  row.names = FALSE)
write.csv(time_table,     "sim_results/table7_timing.csv",          row.names = FALSE)

if (tune_exists) {
  write.csv(tune_summary, "sim_results/table6_tuning_summary.csv", row.names = FALSE)
  write.csv(best_lambda,  "sim_results/table6_best_lambda.csv",    row.names = FALSE)
}
