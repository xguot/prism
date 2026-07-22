#!/usr/bin/env Rscript
# lagrange manuscript figures — exact mda plotting style.
#
# Theme, color, line-type, facet, and sizing conventions are taken verbatim
# from appendix_plots_from_tex.R and plot_results.R.  Every figure uses the
# same theme_mda() and scale_method_aes blocks so the output is visually
# uniform and matches the mda canon.
#
# Output (22 plots → figs/, 7 imputation methods, FIML excluded — it estimates
# parameters but does not produce completed data):
#   fig1a_frobenius_overview_mar.pdf         facet_grid(dist ~ N_label)
#   fig1b_frobenius_overview_mnar.pdf
#   fig1c_frobenius_perdist_mar.pdf          facet_wrap(~ N_label, ncol=3) per dist
#   fig1d_frobenius_perdist_mnar.pdf
#   fig2a_slope_bias_overview_mar.pdf        facet_grid(dist ~ N_label)
#   fig2b_slope_bias_overview_mnar.pdf
#   fig2c_slope_bias_perdist_mar.pdf         facet_wrap(~ N_label, ncol=3) per dist
#   fig2d_slope_bias_perdist_mnar.pdf
#   fig3_relbias_heatmap.pdf                tile heatmap, gradient2, MAR pooled
#   fig4_outlier_degradation.pdf            bar chart, Δ Frobenius Normal→Outlier

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ══════════════════════════════════════════════════════════════════════════════
# Shared theme — exact copy of mda theme
# ══════════════════════════════════════════════════════════════════════════════
theme_mda <- function(legend_pos = "bottom") {
  theme_bw() +
    theme(
      legend.position      = legend_pos,
      legend.text          = element_text(size = 8),
      legend.title         = element_blank(),
      legend.background    = element_blank(),
      legend.key.width     = unit(0.8, "cm"),
      legend.spacing.x     = unit(0.2, "cm"),
      legend.margin        = margin(t = 0, r = 0, b = 0, l = 0),
      strip.text           = element_text(size = 9, face = "bold"),
      axis.title           = element_text(size = 10, face = "bold"),
      axis.text            = element_text(size = 8),
      panel.grid.minor     = element_blank(),
      legend.direction     = "horizontal"
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Shared scales — mda color palette for baselines, warm tones for Smriti
# ══════════════════════════════════════════════════════════════════════════════
method_levels <- c(
  "FIML_Predict", "MICE", "missForest", "missRanger",
  "Smriti_FIML", "Smriti_Default", "Smriti_Robust"
)

method_colors <- c(
  "FIML_Predict"   = "#00BFC4",   # teal — naive FIML completed
  "MICE"           = "#4DAF4A",   # green — multiple imputation
  "missForest"     = "#984EA3",   # purple — ML single imputation
  "missRanger"     = "#FF7F00",   # orange — fast ML baseline
  "Smriti_FIML"    = "#E41A1C",   # red — primary proposed method
  "Smriti_Default" = "#A65628",   # brown — pairwise-target variant
  "Smriti_Robust"  = "#F781BF"    # pink — robust variant
)

# Solid for proposed methods, dashed/dotted for baselines
method_linetypes <- c(
  "FIML_Predict"   = "dotted",
  "MICE"           = "dotdash",
  "missForest"     = "longdash",
  "missRanger"     = "dotted",
  "Smriti_FIML"    = "solid",
  "Smriti_Default" = "solid",
  "Smriti_Robust"  = "solid"
)

scale_method_aes <- list(
  scale_color_manual(values = method_colors, breaks = method_levels),
  scale_linetype_manual(values = method_linetypes, breaks = method_levels),
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
)

# ══════════════════════════════════════════════════════════════════════════════
# Labels and constants
# ══════════════════════════════════════════════════════════════════════════════
dist_levels  <- c("Normal", "t5", "Outlier", "Lognormal")
dist_labels  <- c("Normal", "Student t(5)", "5% Outliers", "Lognormal")

miss_breaks  <- c(0.05, 0.10, 0.15, 0.30)
miss_labels  <- c("5%", "10%", "15%", "30%")

N_levels     <- c("N = 100", "N = 200", "N = 500",
                  "N = 1k", "N = 5k", "N = 10k")

beta_true    <- c(psi_L = 1, psi_S = 1, psi_LS = 0, beta_L = 6, beta_S = 2)
param_names  <- names(beta_true)
param_labels_tex <- c(
  psi_L  = "sigma[L]^2",
  psi_S  = "sigma[S]^2",
  psi_LS = "sigma[LS]",
  beta_L = "beta[L]",
  beta_S = "beta[S]"
)

dir.create("figs", showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Load and preprocess
# ══════════════════════════════════════════════════════════════════════════════
cat("Loading prod_results.rds ...\n")
prod <- readRDS("sim_results/prod_results.rds")
prod <- prod %>% filter(method != "FIML")

prod$method   <- factor(prod$method, levels = method_levels)
prod$dist     <- factor(prod$dist, levels = dist_levels, labels = dist_labels)
prod$N_label  <- factor(
  paste0("N = ", ifelse(prod$N >= 1000, paste0(prod$N / 1000, "k"), prod$N)),
  levels = N_levels
)
prod$miss_pct <- factor(paste0(prod$miss * 100, "%"), levels = miss_labels)

# Aggregate across Monte Carlo replicates
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

# Per-parameter long table (for heatmap)
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
# ── Helper: overview plot (facet_grid) ───────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
make_overview <- function(data, y_var, y_label, title_prefix) {
  p <- ggplot(data, aes(x = miss, y = .data[[y_var]], group = method)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = y_label,
         title = title_prefix) +
    theme_mda()
  p
}

# ── Helper: per-distribution plot (facet_wrap) ───────────────────────────────
make_perdist <- function(data, y_var, y_label, title_template) {
  p <- ggplot(data, aes(x = miss, y = .data[[y_var]], group = method)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_wrap(~ N_label, ncol = 3, scales = "free_y") +
    labs(x = "Missingness Rate", y = y_label,
         title = title_template) +
    theme_mda()
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1a / 1b — Frobenius Distance Overview (MAR, MNAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nFigure 1a: Frobenius Distance Overview — MAR\n")
fig1a <- agg %>% filter(mech == "MAR")
p1a <- make_overview(fig1a, "f_dist_mean",
                     "Frobenius Distance to True Covariance",
                     "Covariance Recovery \u2014 MAR")
ggsave("figs/fig1a_frobenius_overview_mar.pdf",
       p1a, width = 30, height = 20, units = "cm")

cat("Figure 1b: Frobenius Distance Overview — MNAR\n")
fig1b <- agg %>% filter(mech == "MNAR")
p1b <- make_overview(fig1b, "f_dist_mean",
                     "Frobenius Distance to True Covariance",
                     "Covariance Recovery \u2014 MNAR")
ggsave("figs/fig1b_frobenius_overview_mnar.pdf",
       p1b, width = 30, height = 20, units = "cm")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1c / 1d — Frobenius Distance Per-Distribution (MAR, MNAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nFigure 1c: Frobenius Distance Per-Distribution — MAR\n")
for (d in dist_labels) {
  fig <- agg %>% filter(mech == "MAR", dist == d)
  p <- make_perdist(fig, "f_dist_mean",
                    "Frobenius Distance to True Covariance",
                    paste0("Covariance Recovery \u2014 MAR \u2014 ", d))
  safe <- gsub("[ ()%]+", "_", d)
  ggsave(sprintf("figs/fig1c_frobenius_mar_%s.pdf", safe),
         p, width = 24, height = 14, units = "cm")
}

cat("Figure 1d: Frobenius Distance Per-Distribution — MNAR\n")
for (d in dist_labels) {
  fig <- agg %>% filter(mech == "MNAR", dist == d)
  p <- make_perdist(fig, "f_dist_mean",
                    "Frobenius Distance to True Covariance",
                    paste0("Covariance Recovery \u2014 MNAR \u2014 ", d))
  safe <- gsub("[ ()%]+", "_", d)
  ggsave(sprintf("figs/fig1d_frobenius_mnar_%s.pdf", safe),
         p, width = 24, height = 14, units = "cm")
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2a / 2b — Slope Variance Bias Overview (MAR, MNAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nFigure 2a: Slope Variance Bias Overview — MAR\n")
fig2a <- agg %>% filter(mech == "MAR")
p2a <- make_overview(fig2a, "s_var_bias_mean",
                     "Slope Variance Relative Bias (%)",
                     "Parameter Recovery \u2014 MAR") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4)
ggsave("figs/fig2a_slope_bias_overview_mar.pdf",
       p2a, width = 30, height = 20, units = "cm")

cat("Figure 2b: Slope Variance Bias Overview — MNAR\n")
fig2b <- agg %>% filter(mech == "MNAR")
p2b <- make_overview(fig2b, "s_var_bias_mean",
                     "Slope Variance Relative Bias (%)",
                     "Parameter Recovery \u2014 MNAR") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4)
ggsave("figs/fig2b_slope_bias_overview_mnar.pdf",
       p2b, width = 30, height = 20, units = "cm")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2c / 2d — Slope Variance Bias Per-Distribution (MAR, MNAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nFigure 2c: Slope Variance Bias Per-Distribution — MAR\n")
for (d in dist_labels) {
  fig <- agg %>% filter(mech == "MAR", dist == d)
  p <- make_perdist(fig, "s_var_bias_mean",
                    "Slope Variance Relative Bias (%)",
                    paste0("Parameter Recovery \u2014 MAR \u2014 ", d)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4)
  safe <- gsub("[ ()%]+", "_", d)
  ggsave(sprintf("figs/fig2c_slope_bias_mar_%s.pdf", safe),
         p, width = 24, height = 14, units = "cm")
}

cat("Figure 2d: Slope Variance Bias Per-Distribution — MNAR\n")
for (d in dist_labels) {
  fig <- agg %>% filter(mech == "MNAR", dist == d)
  p <- make_perdist(fig, "s_var_bias_mean",
                    "Slope Variance Relative Bias (%)",
                    paste0("Parameter Recovery \u2014 MNAR \u2014 ", d)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4)
  safe <- gsub("[ ()%]+", "_", d)
  ggsave(sprintf("figs/fig2d_slope_bias_mnar_%s.pdf", safe),
         p, width = 24, height = 14, units = "cm")
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Relative Bias Heatmap (gradient2 tile, MAR pooled)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nFigure 3: Relative Bias Heatmap (MAR, pooled)\n")

heatmap_data <- all_params %>%
  filter(mech == "MAR") %>%
  group_by(dist, method, param) %>%
  summarise(RelBias = mean(relbias, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    param_label = factor(param, levels = param_names,
                         labels = param_labels_tex)
  )

p3 <- ggplot(heatmap_data,
       aes(x = param_label, y = method, fill = RelBias)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~ dist, nrow = 1) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 0, name = "Rel Bias (%)") +
  scale_x_discrete(labels = scales::parse_format()) +
  labs(title = "Relative Bias of GCM Parameters by Method and Distribution",
       subtitle = "MAR \u2014 pooled across sample sizes and missingness rates",
       x = "Parameter", y = "Method") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    panel.grid     = element_blank(),
    strip.text     = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("figs/fig3_relbias_heatmap.pdf",
       p3, width = 32, height = 14, units = "cm")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Outlier Degradation (Δ Frobenius Normal → Outlier, MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nFigure 4: Outlier Degradation (MAR)\n")

fig4 <- agg %>%
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
  mutate(miss_pct = factor(paste0(miss * 100, "%"), levels = miss_labels))

p4 <- ggplot(fig4, aes(x = method, y = Delta_mean, fill = method)) +
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

ggsave("figs/fig4_outlier_degradation.pdf",
       p4, width = 27, height = 10, units = "cm")

cat("\nAll figures saved to figs/\n")
