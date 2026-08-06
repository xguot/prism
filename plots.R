suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

prod <- readRDS("sim_results/prod_results.rds")

# Shared theme — adapted from mda/plot_results.R
theme_mda <- function(legend_pos = "bottom") {
  theme_bw() +
    theme(
      legend.position      = legend_pos,
      legend.text          = element_text(size = 7),
      legend.background    = element_blank(),
      legend.key.width     = unit(0.8, "cm"),
      legend.spacing.x     = unit(0.2, "cm"),
      legend.margin        = margin(t = 0, r = 0, b = 0, l = 0),
      strip.text           = element_text(size = 8, face = "bold"),
      axis.title           = element_text(size = 10),
      axis.text            = element_text(size = 7),
      panel.grid.minor     = element_blank(),
      legend.direction     = "horizontal"
    )
}

# Method levels and aesthetics
method_levels <- c("FIML", "FIML_lavPredict", "FIML_stochastic",
                   "MICE_norm", "missForest", "missRanger",
                   "PRISM", "PRISM_MI")
method_colors <- c(
  "FIML"             = "#000000", "FIML_lavPredict"  = "#999999",
  "FIML_stochastic"  = "#56B4E9", "MICE_norm"        = "#E69F00",
  "missForest"       = "#009E73", "missRanger"       = "#0072B2",
  "PRISM"            = "#D55E00", "PRISM_MI"         = "#CC79A7"
)
method_shapes <- c(
  "FIML" = 17, "FIML_lavPredict" = 15, "FIML_stochastic" = 18,
  "MICE_norm" = 16, "missForest" = 8, "missRanger" = 4,
  "PRISM" = 19, "PRISM_MI" = 1
)
scale_method_aes <- list(
  scale_color_manual(values = method_colors, breaks = method_levels),
  scale_shape_manual(values = method_shapes, breaks = method_levels),
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2)),
         shape = guide_legend(nrow = 2))
)

# Factor levels
dist_levels <- c("Normal", "t5", "Outlier", "Lognormal")
miss_breaks <- c(0.05, 0.10, 0.15, 0.30)
miss_labels <- c("5%", "10%", "15%", "30%")

prod$method <- factor(prod$method, levels = method_levels)
prod$dist   <- factor(prod$dist, levels = dist_levels)
prod$N_label <- factor(paste0("N = ", prod$N),
                       levels = paste0("N = ", c(100, 200, 500, 1000, 5000, 10000)))

# Aggregate Frobenius
frob_agg <- prod %>%
  filter(is.finite(f_dist)) %>%
  group_by(N, N_label, miss, dist, mech, method) %>%
  summarise(frob_mean = mean(f_dist, na.rm = TRUE),
            frob_sd   = sd(f_dist, na.rm = TRUE), .groups = "drop")

# Per-parameter bias long table
beta_true   <- c(psi_L = 1, psi_S = 1, psi_LS = 0, beta_L = 6, beta_S = 2)
param_names <- names(beta_true)
param_labels <- c(psi_L  = "Var Intercept", psi_S = "Var Slope",
                  psi_LS = "Covariance", beta_L = "Intercept", beta_S = "Slope")
param_files <- c(psi_L  = "Var_Intercept", psi_S = "Var_Slope",
                 psi_LS = "Covariance", beta_L = "Intercept", beta_S = "Slope")
est_cols <- c(psi_L = "est_var_L", psi_S = "est_var_S", psi_LS = "est_cov_LS",
              beta_L = "est_L", beta_S = "est_S")

bias_agg <- do.call(rbind, lapply(param_names, function(pn) {
  tv <- beta_true[pn]; ec <- est_cols[pn]
  prod %>%
    filter(is.finite(.data[[ec]])) %>%
    mutate(param = pn, est = .data[[ec]],
           rel_bias = if (abs(tv) < 1e-12) est - tv else 100 * (est - tv) / tv) %>%
    select(N, N_label, miss, dist, mech, method, param, rel_bias)
})) %>%
  group_by(N, N_label, miss, dist, mech, method, param) %>%
  summarise(arb_mean = mean(abs(rel_bias), na.rm = TRUE),
            arb_sd   = sd(abs(rel_bias), na.rm = TRUE), .groups = "drop")

dir.create("figs", showWarnings = FALSE)

# FIGURE 1 — Frobenius MAR
cat("Plotting Frobenius (MAR)...\n")
p <- frob_agg %>% filter(mech == "MAR") %>%
  ggplot(aes(x = miss, y = frob_mean, group = method)) +
  geom_line(aes(color = method), linewidth = 0.5) +
  geom_point(aes(color = method, shape = method), size = 1.0) +
  scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
  scale_method_aes +
  facet_grid(dist ~ N_label, scales = "free_y") +
  labs(x = "Missingness Rate", y = "Frobenius Distance", title = "Covariance Recovery: Frobenius Distance to True Covariance (MAR)") +
  theme_mda()
ggsave("figs/Frobenius_Distance_MAR.pdf", p, width = 30, height = 20, units = "cm")

# FIGURE 2 — Frobenius MNAR
cat("Plotting Frobenius (MNAR)...\n")
p <- frob_agg %>% filter(mech == "MNAR") %>%
  ggplot(aes(x = miss, y = frob_mean, group = method)) +
  geom_line(aes(color = method), linewidth = 0.5) +
  geom_point(aes(color = method, shape = method), size = 1.0) +
  scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
  scale_method_aes +
  facet_grid(dist ~ N_label, scales = "free_y") +
  labs(x = "Missingness Rate", y = "Frobenius Distance", title = "Covariance Recovery: Frobenius Distance to True Covariance (MNAR)") +
  theme_mda()
ggsave("figs/Frobenius_Distance_MNAR.pdf", p, width = 30, height = 20, units = "cm")

# FIGURES 3-7 — Absolute Relative Bias per parameter (MAR)
for (pn in param_names) {
  cat(sprintf("Plotting Bias: %s...\n", pn))
  tv <- beta_true[pn]
  y_lab <- if (abs(tv) < 1e-12) "Absolute Raw Bias" else paste0("Absolute Relative Bias in ", param_labels[pn], " (%)")
  hline <- if (abs(tv) < 1e-12) NULL else
    geom_hline(yintercept = 10, linetype = "dashed", color = "grey50", linewidth = 0.4)

  fig <- bias_agg %>% filter(mech == "MAR", param == pn)
  p <- ggplot(fig, aes(x = miss, y = arb_mean, group = method)) +
    hline +
    geom_line(aes(color = method), linewidth = 0.5) +
    geom_point(aes(color = method, shape = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = y_lab,
         title = param_labels[pn]) +
    theme_mda()
  ggsave(sprintf("figs/bias_mar_%s.pdf", param_files[pn]), p, width = 30, height = 20, units = "cm")
}

cat(sprintf("\nDone. %d figures saved to figs/\n", 2 + length(param_names)))
