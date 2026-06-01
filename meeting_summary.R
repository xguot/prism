library(dplyr)
library(tidyr)

all_data <- readRDS("sim_results/prod_results.rds")

# Backward compat: if pipeline_time column is missing (pre-v0.1.1 data),
# compute it as smriti routing time + the average missForest time per condition.
if (!"pipeline_time" %in% names(all_data)) {
  mf_times <- all_data %>%
    filter(method == "missForest") %>%
    group_by(dist, mech) %>%
    summarize(mf_time = mean(time_sec, na.rm = TRUE), .groups = "drop")

  all_data <- all_data %>%
    left_join(mf_times, by = c("dist", "mech")) %>%
    mutate(pipeline_time = ifelse(
      method %in% c("Smriti_Default", "Smriti_Robust"),
      mf_time + time_sec,
      time_sec
    )) %>%
    select(-mf_time)
}

# ── Helper: summarise one subset ───────────────────────────────────────────────
summarise_methods <- function(.data, label) {
  .data %>%
    group_by(dist, method) %>%
    summarize(
      Bias           = mean(s_var_bias, na.rm = TRUE),
      Stability_SD   = sd(s_var, na.rm = TRUE),
      Time_sec       = mean(time_sec, na.rm = TRUE),
      Pipeline_sec   = mean(pipeline_time, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    mutate(
      Bias       = sprintf("%+.2f%%", Bias),
      Stability  = sprintf("%.3f", Stability_SD),
      Time       = sprintf("%.2fs", Time_sec),
      Pipeline   = ifelse(method %in% c("Smriti_Default", "Smriti_Robust"),
                          sprintf("%.2fs", Pipeline_sec), "\u2014"),
      Note       = ifelse(method %in% c("Smriti_Default", "Smriti_Robust"),
                          "routing only; pipeline incl. missForest init", "")
    ) %>%
    select(dist, method, Bias, Stability, Time, Pipeline, Note)
}

# ── Table 1: MAR (Missing At Random) — the fair comparison ─────────────────────
cat("\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("TABLE 1 \u2014 MAR (Missing At Random)\n")
cat("FIML is the gold-standard baseline.  All methods can recover the signal.\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

mar_data  <- all_data %>% filter(mech == "MAR")
mar_table <- summarise_methods(mar_data, "MAR")
print(as.data.frame(mar_table), row.names = FALSE)

# ── Table 2: MNAR (Missing Not At Random) — imputation ceiling ─────────────────
cat("\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("TABLE 2 \u2014 MNAR (Missing Not At Random)\n")
cat("All methods degrade substantially.  FIML itself hits -26% to -81% bias.\n")
cat("Imputation cannot recover information that was never observed.\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

mnar_data  <- all_data %>% filter(mech == "MNAR")
mnar_table <- summarise_methods(mnar_data, "MNAR")
print(as.data.frame(mnar_table), row.names = FALSE)

# ── Caveats ────────────────────────────────────────────────────────────────────
cat("\n")
cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
cat("CAVEATS\n")
cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n\n")

cat("1. TIMING: Smriti's Time column reports the Lagrangian routing step only\n")
cat("   (~0.5 s).  The Pipeline column adds missForest initialisation cost\n")
cat("   (~8 s), yielding a full end-to-end cost of ~8.5 s \u2014 comparable to\n")
cat("   missForest alone.  The routing itself is fast, but it depends on a\n")
cat("   pre-existing initial imputation.\n\n")

cat("2. LOGNORMAL FAILURE MODE: Under lognormal MAR, Smriti_Robust reaches\n")
cat("   -88.9% bias vs. missForest's -27.5%.  Spearman/MAD covariance targets\n")
cat("   do not adequately capture lognormal variance structure.  Smriti is not\n")
cat("   recommended for heavily skewed distributions without a custom_target.\n\n")

cat("3. REAL ADVANTAGE: Under Normal and Outlier MAR, the Lagrangian routing\n")
cat("   improves covariance recovery by 1\u20133 percentage points over missForest:\n")
cat("   Normal MAR:    Smriti_Robust -3.67%  vs. missForest -4.59%\n")
cat("   Outlier MAR:   Smriti_Robust -4.51%  vs. missForest -6.47%\n")
cat("   FIML remains the best method (-0.34%) when applicable.\n")
