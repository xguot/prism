# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Plots for smriti
# Style adapted from mda/analysis.R — heatmaps, MSE bars, cleaner facets
# Uses prod_results.rds
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Shared plot theme ─────────────────────────────────────────────────────────
theme_smriti <- function(legend_pos = "bottom") {
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

# ── Shared colour / linetype scale ───────────────────────────────────────────
method_levels <- c("FIML", "FIML_Predict", "MICE", "missForest", "missRanger",
                   "Smriti_FIML", "Smriti_Default", "Smriti_Robust")
method_colors <- c(
  "FIML"           = "#999999",
  "FIML_Predict"   = "#B0B0B0",
  "MICE"           = "#56B4E9",
  "missForest"     = "#009E73",
  "missRanger"     = "#0072B5",
  "Smriti_FIML"    = "#E6A000",
  "Smriti_Default" = "#D55E00",
  "Smriti_Robust"  = "#CC0000"
)
method_linetypes <- c(
  "FIML"           = "dotted",
  "FIML_Predict"   = "dotted",
  "MICE"           = "dashed",
  "missForest"     = "dashed",
  "missRanger"     = "dotted",
  "Smriti_FIML"    = "solid",
  "Smriti_Default" = "solid",
  "Smriti_Robust"  = "solid"
)

scale_method_aes <- list(
  scale_color_manual(values = method_colors, breaks = method_levels),
  scale_linetype_manual(values = method_linetypes, breaks = method_levels)
)

# ── Labels ────────────────────────────────────────────────────────────────────
dist_levels  <- c("Normal", "t5", "Outlier", "Lognormal")
dist_labels  <- c("Normal", "Student t(5)", "5% Outliers", "Lognormal")
miss_levels  <- c("5%", "10%", "15%", "30%")
N_levels     <- c("N = 100", "N = 200", "N = 500", "N = 1k", "N = 5k", "N = 10k")

beta_true   <- c(psi_L = 1, psi_S = 1, psi_LS = 0, beta_L = 6, beta_S = 2)
param_names <- names(beta_true)
param_labels_tex <- c(
  psi_L  = "sigma[L]^2",
  psi_S  = "sigma[S]^2",
  psi_LS = "sigma[LS]",
  beta_L = "beta[L]",
  beta_S = "beta[S]"
)

dir.create("manuscript_figures", showWarnings = FALSE, recursive = TRUE)

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

# Aggregate
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

# ── Per-parameter long table for heatmap / MSE plots ────────────────────────
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
    select(N, miss, dist, mech, method, param, est, bias_raw, relbias)
}))

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Frobenius Distance by Missingness Rate (MAR, headline line plot)
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
# FIGURE 2 — Relative Bias Heatmap (mda-style, MAR, pooled across conditions)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 2: Relative Bias Heatmap (MAR, pooled)\n")

heatmap_data <- all_params %>%
  filter(mech == "MAR") %>%
  group_by(dist, method, param) %>%
  summarise(RelBias = mean(relbias, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    param_label = factor(param, levels = param_names,
                         labels = param_labels_tex)
  )

p2 <- ggplot(heatmap_data,
       aes(x = param_label, y = method, fill = RelBias)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ dist, nrow = 1) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 0, name = "Rel Bias (%)") +
  scale_x_discrete(labels = scales::parse_format()) +
  labs(title = "Relative Bias of GCM Parameters by Method and Distribution",
       subtitle = "MAR, pooled across sample sizes and missingness rates",
       x = "Parameter", y = "Method") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    panel.grid     = element_blank(),
    strip.text     = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("manuscript_figures/fig2_relbias_heatmap.pdf",
       plot = p2, width = 32, height = 14, units = "cm")
ggsave("manuscript_figures/fig2_relbias_heatmap.png",
       plot = p2, width = 32, height = 14, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Outlier Degradation (Δ Frobenius: Normal -> Outlier, MAR)
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
# FIGURE 4 — MSE Bar Chart (mda-style, MAR, pooled, faceted by parameter)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 4: MSE Bar Chart by Parameter and Distribution (MAR)\n")

mse_data <- all_params %>%
  filter(mech == "MAR") %>%
  group_by(dist, method, param) %>%
  summarise(
    bias_raw = mean(est, na.rm = TRUE) - beta_true[param[1]],
    ESE      = sd(est, na.rm = TRUE),
    MSE      = bias_raw^2 + ESE^2,
    .groups  = "drop"
  ) %>%
  mutate(
    param_label = factor(param, levels = param_names,
                         labels = param_labels_tex)
  )

p4 <- ggplot(mse_data, aes(x = method, y = MSE, fill = method)) +
  geom_col() +
  scale_fill_manual(values = method_colors, guide = "none") +
  facet_grid(param_label ~ dist, scales = "free_y",
             labeller = labeller(param_label = label_parsed)) +
  labs(title = "MSE of GCM Parameters by Method and Distribution",
       subtitle = "MAR, pooled across sample sizes and missingness rates",
       x = "", y = "MSE") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    strip.text.y   = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave("manuscript_figures/fig4_mse_bars.pdf",
       plot = p4, width = 28, height = 20, units = "cm")
ggsave("manuscript_figures/fig4_mse_bars.png",
       plot = p4, width = 28, height = 20, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Slope Variance Bias by Missingness Rate (MAR + MNAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 5: Slope Variance Bias by Missingness (MAR & MNAR)\n")

fig5 <- agg %>% filter(mech %in% c("MAR", "MNAR"))

p5 <- ggplot(fig5, aes(x = miss, y = s_var_bias_mean, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_continuous(
    breaks = c(0.05, 0.10, 0.15, 0.30),
    labels = c("5%", "10%", "15%", "30%")
  ) +
  scale_method_aes +
  facet_grid(dist ~ mech + N_label) +
  labs(x = "Missingness Rate",
       y = "Slope Variance Relative Bias (%)") +
  theme_smriti()

ggsave("manuscript_figures/fig5_slope_bias.pdf",
       plot = p5, width = 52, height = 18, units = "cm")
ggsave("manuscript_figures/fig5_slope_bias.png",
       plot = p5, width = 52, height = 18, units = "cm", dpi = 300)

cat("\nAll figures saved to manuscript_figures/\n")
