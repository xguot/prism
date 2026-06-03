# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Plots for smriti
# Uses existing prod_results.rds — no re-simulation needed.
# Plot conventions adapted from plot_demo (Tang & Tong, UVA) simulation study:
#   theme_bw + facet_grid, colour/linetype per method, 27×12-15 cm output.
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Shared plot theme (plot_demo convention) ──────────────────────────────────
theme_smriti <- function(legend_pos = c(0.85, 0.75)) {
  theme_bw() +
    theme(
      legend.position      = legend_pos,
      legend.text          = element_text(size = 8),
      legend.background    = element_blank(),
      legend.key.width     = unit(1.2, "cm"),
      strip.text           = element_text(size = 9, face = "bold"),
      axis.title           = element_text(size = 10),
      axis.text            = element_text(size = 8),
      panel.grid.minor     = element_blank()
    )
}

# ── Shared colour / linetype scale (plot_demo convention) ─────────────────────
method_levels <- c("FIML", "MICE", "missForest", "missRanger",
                   "Smriti_Default", "Smriti_Robust")
method_colors <- c(
  "FIML"           = "#999999",
  "MICE"           = "#56B4E9",
  "missForest"     = "#009E73",
  "missRanger"     = "#0072B5",
  "Smriti_Default" = "#D55E00",
  "Smriti_Robust"  = "#CC0000"
)
method_linetypes <- c(
  "FIML"           = "dotted",
  "MICE"           = "dashed",
  "missForest"     = "dashed",
  "missRanger"     = "dotted",
  "Smriti_Default" = "solid",
  "Smriti_Robust"  = "solid"
)

scale_method_aes <- list(
  scale_color_manual(values = method_colors, breaks = method_levels),
  scale_linetype_manual(values = method_linetypes, breaks = method_levels)
)

dist_levels  <- c("Normal", "t5", "Outlier", "Lognormal")
dist_labels  <- c("Normal", "Student t(5)", "5% Outliers", "Lognormal")
miss_levels  <- c("5%", "10%", "15%", "30%")
N_levels     <- c("N = 100", "N = 200", "N = 500", "N = 1k", "N = 5k", "N = 10k")

# ── Parameter truths for reference lines ──────────────────────────────────────
param_truths <- c(
  est_L     = 6,  est_S     = 2,
  est_var_L = 1,  est_var_S = 1,  est_cov_LS = 0,
  bias_L    = 0,  bias_S    = 0,
  bias_var_L = 0, bias_var_S = 0, bias_cov_LS = 0
)

# ══════════════════════════════════════════════════════════════════════════════
# Load & preprocess
# ══════════════════════════════════════════════════════════════════════════════
prod <- readRDS("sim_results/prod_results.rds")
prod$method   <- factor(prod$method, levels = method_levels)
prod$dist     <- factor(prod$dist, levels = dist_levels, labels = dist_labels)
prod$N_label  <- factor(paste0("N = ", ifelse(prod$N >= 1000,
                           paste0(prod$N / 1000, "k"), prod$N)),
                        levels = N_levels)
prod$miss_pct <- factor(paste0(prod$miss * 100, "%"), levels = miss_levels)

# ── Aggregate ─────────────────────────────────────────────────────────────────
agg <- prod %>%
  group_by(N, N_label, miss, miss_pct, dist, mech, method) %>%
  summarise(
    across(c(f_dist, s_var, s_var_bias, s_se,
             est_L, est_S, est_var_L, est_var_S, est_cov_LS,
             bias_L, bias_S, bias_var_L, bias_var_S, bias_cov_LS),
           list(mean = ~ mean(.x, na.rm = TRUE),
                sd   = ~ sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    time_mean = mean(time_sec, na.rm = TRUE),
    .groups   = "drop"
  )

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Frobenius Distance by Missingness Rate (MAR, headline)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 1: Frobenius Distance by Missingness Rate (MAR)\n")

fig1 <- agg %>% filter(mech == "MAR")

p1 <- ggplot(fig1, aes(x = miss, y = f_dist_mean, group = method)) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_continuous(
    breaks = c(0.05, 0.10, 0.15, 0.30),
    labels = c("5%", "10%", "15%", "30%")
  ) +
  scale_method_aes +
  facet_grid(dist ~ N_label) +
  labs(x = "Missingness Rate",
       y = "Frobenius Distance to True Covariance") +
  theme_smriti()

ggsave("manuscript_figures/fig1_frobenius_by_miss.pdf",
       plot = p1, width = 27, height = 18, units = "cm")
ggsave("manuscript_figures/fig1_frobenius_by_miss.png",
       plot = p1, width = 27, height = 18, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Frobenius Distance by Sample Size (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 2: Frobenius Distance by Sample Size (MAR)\n")

p2 <- ggplot(fig1, aes(x = N, y = f_dist_mean, group = method)) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000, 10000),
    labels = c("100", "200", "500", "1k", "5k", "10k")
  ) +
  scale_method_aes +
  facet_grid(dist ~ miss_pct) +
  labs(x = "Sample Size (N)",
       y = "Frobenius Distance to True Covariance") +
  theme_smriti()

ggsave("manuscript_figures/fig2_frobenius_by_N.pdf",
       plot = p2, width = 27, height = 18, units = "cm")
ggsave("manuscript_figures/fig2_frobenius_by_N.png",
       plot = p2, width = 27, height = 18, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Outlier Degradation (Δ Frobenius: Normal → Outlier, MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 3: Outlier Degradation (MAR)\n")

fig3 <- agg %>%
  filter(mech == "MAR", dist %in% c("Normal", "5% Outliers")) %>%
  select(N, miss, dist, method, f_dist_mean) %>%
  pivot_wider(names_from = dist, values_from = f_dist_mean) %>%
  rename(Normal = "Normal", Outlier = "5% Outliers") %>%
  mutate(Delta = Outlier - Normal) %>%
  group_by(miss, method) %>%
  summarise(
    Delta_mean = mean(Delta, na.rm = TRUE),
    Delta_sd   = sd(Delta, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(miss_pct = factor(paste0(miss * 100, "%"), levels = miss_levels))

p3 <- ggplot(fig3, aes(x = method, y = Delta_mean, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Delta_mean - Delta_sd, ymax = Delta_mean + Delta_sd),
                width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = method_colors, guide = "none") +
  facet_wrap(~ miss_pct, nrow = 1) +
  labs(x = "",
       y = expression(Delta ~ "Frobenius (Outlier - Normal)")) +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    strip.text         = element_text(size = 9, face = "bold"),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )

ggsave("manuscript_figures/fig3_outlier_degradation.pdf",
       plot = p3, width = 27, height = 10, units = "cm")
ggsave("manuscript_figures/fig3_outlier_degradation.png",
       plot = p3, width = 27, height = 10, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Parameter Bias: Relative Bias by N (MAR, plot_demo style)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 4: Parameter Relative Bias (MAR)\n")

bias_long <- agg %>%
  filter(mech == "MAR", miss == 0.15) %>%
  select(N, dist, method, bias_L_mean, bias_S_mean,
         bias_var_L_mean, bias_var_S_mean, bias_cov_LS_mean) %>%
  pivot_longer(
    cols      = c(bias_L_mean, bias_S_mean, bias_var_L_mean, bias_var_S_mean, bias_cov_LS_mean),
    names_to  = "Parameter",
    values_to = "Bias"
  ) %>%
  mutate(Parameter = factor(Parameter,
    levels = c("bias_L_mean", "bias_S_mean", "bias_var_L_mean", "bias_var_S_mean", "bias_cov_LS_mean"),
    labels = c("beta[L]", "beta[S]", "sigma[L]^2", "sigma[S]^2", "sigma[LS]")
  ))

p4 <- ggplot(bias_long, aes(x = N, y = abs(Bias), group = method)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.2) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000, 10000),
    labels = c("100", "200", "500", "1k", "5k", "10k")
  ) +
  scale_method_aes +
  facet_grid(dist ~ Parameter, labeller = label_parsed) +
  labs(x = "Sample Size (N)",
       y = "|Relative Bias| (%)",
       subtitle = "Dashed line: 10% acceptable bias threshold.  15% missingness, MAR.") +
  theme_smriti()

ggsave("manuscript_figures/fig4_parameter_bias.pdf",
       plot = p4, width = 30, height = 18, units = "cm")
ggsave("manuscript_figures/fig4_parameter_bias.png",
       plot = p4, width = 30, height = 18, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Slope Variance Bias (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 5: Slope Variance Bias (MAR)\n")

fig5 <- agg %>% filter(mech == "MAR")

p5 <- ggplot(fig5, aes(x = N, y = s_var_bias_mean, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000, 10000),
    labels = c("100", "200", "500", "1k", "5k", "10k")
  ) +
  scale_method_aes +
  facet_grid(dist ~ miss_pct) +
  labs(x = "Sample Size (N)",
       y = expression("Bias in Recovered " * sigma[s]^2 * " (%)")) +
  theme_smriti()

ggsave("manuscript_figures/fig5_slope_bias.pdf",
       plot = p5, width = 27, height = 18, units = "cm")
ggsave("manuscript_figures/fig5_slope_bias.png",
       plot = p5, width = 27, height = 18, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Timing Comparison (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 6: Timing Comparison (MAR)\n")

fig6 <- agg %>%
  filter(mech == "MAR") %>%
  group_by(method) %>%
  summarise(
    Time_mean = mean(time_mean, na.rm = TRUE),
    Time_sd   = sd(time_mean, na.rm = TRUE),
    .groups   = "drop"
  )

p6 <- ggplot(fig6, aes(x = reorder(method, Time_mean), y = Time_mean, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, Time_mean - Time_sd), ymax = Time_mean + Time_sd),
                width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(x = "", y = "Mean Time (seconds)",
       subtitle = "Error bars: ± 1 SD.  Smriti routing step only (~0.6 s).") +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor   = element_blank()
  )

ggsave("manuscript_figures/fig6_timing.pdf",
       plot = p6, width = 18, height = 10, units = "cm")
ggsave("manuscript_figures/fig6_timing.png",
       plot = p6, width = 18, height = 10, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — t5 vs Normal: Heavy-tailed Robustness (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 7: Heavy-tailed (t5) Comparison (MAR)\n")

fig7 <- agg %>%
  filter(mech == "MAR", dist %in% c("Normal", "Student t(5)"))

p7 <- ggplot(fig7, aes(x = N, y = f_dist_mean, group = method)) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000, 10000),
    labels = c("100", "200", "500", "1k", "5k", "10k")
  ) +
  scale_method_aes +
  facet_grid(dist ~ miss_pct) +
  labs(x = "Sample Size (N)",
       y = "Frobenius Distance to True Covariance") +
  theme_smriti()

ggsave("manuscript_figures/fig7_t5_comparison.pdf",
       plot = p7, width = 27, height = 12, units = "cm")
ggsave("manuscript_figures/fig7_t5_comparison.png",
       plot = p7, width = 27, height = 12, units = "cm", dpi = 300)

cat("\nAll 7 figures saved to manuscript_figures/\n")
