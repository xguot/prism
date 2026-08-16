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

# Relabel MICE to reflect the rf (ranger) backend
prod$method[prod$method == "MICE"] <- "MICE (RF)"

# Method levels and aesthetics
method_levels <- c("FIML", "MICE (RF)", "missForest", "PRISM", "PRISM_MI")
method_colors <- c(
  "FIML"             = "#000000",
  "MICE (RF)"        = "#E69F00",
  "missForest"       = "#009E73",
  "PRISM"            = "#D55E00", "PRISM_MI"         = "#CC79A7"
)
method_shapes <- c(
  "FIML" = 17, "MICE (RF)" = 16,
  "missForest" = 8, "PRISM" = 19, "PRISM_MI" = 1
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

# Drop methods kept only in supplementary material
keep_methods <- c("FIML", "MICE (RF)", "missForest", "PRISM", "PRISM_MI")

# Aggregate Frobenius
frob_agg <- prod %>%
  filter(is.finite(f_dist), method %in% keep_methods) %>%
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
se_cols <- c(psi_L = "se_var_L", psi_S = "se_var_S", psi_LS = "se_cov_LS",
             beta_L = "se_L", beta_S = "se_S")

bias_agg <- do.call(rbind, lapply(param_names, function(pn) {
  tv <- beta_true[pn]; ec <- est_cols[pn]
  prod %>%
    filter(is.finite(.data[[ec]]), method %in% keep_methods) %>%
    mutate(param = pn, est = .data[[ec]],
           rel_bias = if (abs(tv) < 1e-12) est - tv else 100 * (est - tv) / tv) %>%
    select(N, N_label, miss, dist, mech, method, param, rel_bias)
})) %>%
  group_by(N, N_label, miss, dist, mech, method, param) %>%
  summarise(arb_mean = mean(abs(rel_bias), na.rm = TRUE),
            arb_sd   = sd(abs(rel_bias), na.rm = TRUE), .groups = "drop")

# Aggregate SE ratio
se_ratio_agg <- do.call(rbind, lapply(param_names, function(pn) {
  sc <- se_cols[pn]
  if (!sc %in% names(prod)) return(NULL)
  prod %>%
    filter(is.finite(.data[[sc]]), method %in% keep_methods) %>%
    group_by(N, N_label, miss, dist, mech, method) %>%
    summarise(
      empirical_se = sd(.data[[est_cols[pn]]], na.rm = TRUE),
      avg_model_se = mean(.data[[sc]], na.rm = TRUE),
      se_ratio     = avg_model_se / empirical_se,
      .groups      = "drop"
    ) %>%
    mutate(param = pn)
}))

# Aggregate convergence rate
conv_agg <- prod %>%
  filter(method %in% keep_methods) %>%
  group_by(N, N_label, miss, dist, mech, method) %>%
  summarise(conv_rate = mean(converged, na.rm = TRUE), .groups = "drop")

# Aggregate timing
timing_agg <- prod %>%
  filter(method %in% keep_methods, is.finite(time_sec)) %>%
  group_by(method) %>%
  summarise(mean_time = mean(time_sec, na.rm = TRUE),
            sd_time   = sd(time_sec, na.rm = TRUE), .groups = "drop")

# Aggregate coverage (approx: est +/- 1.96 * se)
cov_agg <- do.call(rbind, lapply(param_names, function(pn) {
  tv <- beta_true[pn]; ec <- est_cols[pn]; sc <- se_cols[pn]
  if (!sc %in% names(prod)) return(NULL)
  prod %>%
    filter(is.finite(.data[[ec]]), is.finite(.data[[sc]]), method %in% keep_methods) %>%
    mutate(
      lo = .data[[ec]] - 1.96 * .data[[sc]],
      hi = .data[[ec]] + 1.96 * .data[[sc]],
      covered = (lo <= tv) & (hi >= tv)
    ) %>%
    group_by(N, N_label, miss, dist, mech, method) %>%
    summarise(coverage = mean(covered, na.rm = TRUE), .groups = "drop") %>%
    mutate(param = pn)
}))

# Aggregate MSE (bias^2 + empirical variance)
mse_agg <- do.call(rbind, lapply(param_names, function(pn) {
  tv <- beta_true[pn]; ec <- est_cols[pn]
  prod %>%
    filter(is.finite(.data[[ec]]), method %in% keep_methods) %>%
    group_by(N, N_label, miss, dist, mech, method) %>%
    summarise(
      mse = (mean(.data[[ec]], na.rm = TRUE) - tv)^2 + var(.data[[ec]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(param = pn)
}))

dir.create("figs", showWarnings = FALSE)

# FIGURE 1 — Frobenius MAR
cat("Plotting Frobenius (MAR)...\n")
p <- frob_agg %>% filter(mech == "MAR") %>%
  ggplot(aes(x = miss, y = frob_mean, group = method)) +
  geom_line(aes(color = method), linewidth = 0.5) +
  geom_point(aes(color = method, shape = method), size = 1.0) +
  scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
  scale_y_log10() +
  scale_method_aes +
  facet_grid(dist ~ N_label, scales = "free_y") +
  labs(x = "Missingness Rate", y = "Frobenius Distance (log10)", title = "Frobenius Distance to True Covariance (MAR)") +
  theme_mda()
ggsave("figs/Frobenius_Distance_MAR.pdf", p, width = 30, height = 20, units = "cm")

# FIGURE 2 — Frobenius MNAR
cat("Plotting Frobenius (MNAR)...\n")
p <- frob_agg %>% filter(mech == "MNAR") %>%
  ggplot(aes(x = miss, y = frob_mean, group = method)) +
  geom_line(aes(color = method), linewidth = 0.5) +
  geom_point(aes(color = method, shape = method), size = 1.0) +
  scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
  scale_y_log10() +
  scale_method_aes +
  facet_grid(dist ~ N_label, scales = "free_y") +
  labs(x = "Missingness Rate", y = "Frobenius Distance (log10)", title = "Frobenius Distance to True Covariance (MNAR)") +
  theme_mda()
ggsave("figs/Frobenius_Distance_MNAR.pdf", p, width = 30, height = 20, units = "cm")

# Absolute Relative Bias per parameter (MAR + MNAR)
for (mech_i in c("MAR", "MNAR")) {
  for (pn in param_names) {
    cat(sprintf("Plotting Bias (%s): %s...\n", mech_i, pn))
    tv <- beta_true[pn]
    y_lab <- if (abs(tv) < 1e-12) "Absolute Raw Bias" else paste0("Absolute Relative Bias in ", param_labels[pn], " (%)")
    hline <- if (abs(tv) < 1e-12) NULL else
      geom_hline(yintercept = 10, linetype = "dashed", color = "grey50", linewidth = 0.4)

    fig <- bias_agg %>% filter(mech == mech_i, param == pn)
    p <- ggplot(fig, aes(x = miss, y = arb_mean, group = method)) +
      hline +
      geom_line(aes(color = method), linewidth = 0.5) +
      geom_point(aes(color = method, shape = method), size = 1.0) +
      scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
      scale_method_aes +
      facet_grid(dist ~ N_label, scales = "free_y") +
      labs(x = "Missingness Rate", y = y_lab,
           title = paste0(param_labels[pn], " (", mech_i, ")")) +
      theme_mda()
    ggsave(sprintf("figs/bias_%s_%s.pdf", tolower(mech_i), param_files[pn]), p,
           width = 30, height = 20, units = "cm")
  }
}

# SE metrics per parameter (MAR + MNAR)
total_figs <- 2 + 4 * length(param_names) * 2

if (!is.null(se_ratio_agg)) {
  for (mech_i in c("MAR", "MNAR")) {
    for (pn in param_names) {
      fig <- se_ratio_agg %>% filter(mech == mech_i, param == pn)

      # SE Ratio
      p <- ggplot(fig, aes(x = miss, y = se_ratio, group = method)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
        geom_line(aes(color = method), linewidth = 0.5) +
        geom_point(aes(color = method, shape = method), size = 1.0) +
        scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
        scale_method_aes +
        facet_grid(dist ~ N_label, scales = "free_y") +
        labs(x = "Missingness Rate", y = "SE Ratio (Model SE / Empirical SE)",
             title = paste0("SE Ratio: ", param_labels[pn], " (", mech_i, ")")) +
        theme_mda()
      ggsave(sprintf("figs/se_ratio_%s_%s.pdf", tolower(mech_i), param_files[pn]), p,
             width = 30, height = 20, units = "cm")

      # Empirical SE
      p <- ggplot(fig, aes(x = miss, y = empirical_se, group = method)) +
        geom_line(aes(color = method), linewidth = 0.5) +
        geom_point(aes(color = method, shape = method), size = 1.0) +
        scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
        scale_y_log10() +
        scale_method_aes +
        facet_grid(dist ~ N_label, scales = "free_y") +
        labs(x = "Missingness Rate", y = "Empirical SE (log10)",
             title = paste0("Empirical SE: ", param_labels[pn], " (", mech_i, ")")) +
        theme_mda()
      ggsave(sprintf("figs/empirical_se_%s_%s.pdf", tolower(mech_i), param_files[pn]), p,
             width = 30, height = 20, units = "cm")

      # Average Model SE
      p <- ggplot(fig, aes(x = miss, y = avg_model_se, group = method)) +
        geom_line(aes(color = method), linewidth = 0.5) +
        geom_point(aes(color = method, shape = method), size = 1.0) +
        scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
        scale_y_log10() +
        scale_method_aes +
        facet_grid(dist ~ N_label, scales = "free_y") +
        labs(x = "Missingness Rate", y = "Average Model SE (log10)",
             title = paste0("Average Model SE: ", param_labels[pn], " (", mech_i, ")")) +
        theme_mda()
      ggsave(sprintf("figs/model_se_%s_%s.pdf", tolower(mech_i), param_files[pn]), p,
             width = 30, height = 20, units = "cm")
    }
  }
}

# Convergence rate (MAR + MNAR)
for (mech_i in c("MAR", "MNAR")) {
  p <- conv_agg %>% filter(mech == mech_i) %>%
    ggplot(aes(x = miss, y = conv_rate, group = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(aes(color = method), linewidth = 0.5) +
    geom_point(aes(color = method, shape = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label) +
    labs(x = "Missingness Rate", y = "Convergence Rate",
         title = paste0("Convergence Rate (", mech_i, ")")) +
    theme_mda()
  ggsave(sprintf("figs/conv_rate_%s.pdf", tolower(mech_i)), p,
         width = 30, height = 20, units = "cm")
}

# Timing comparison (bar)
p <- timing_agg %>%
  ggplot(aes(x = method, y = mean_time, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = pmax(mean_time - sd_time, 1e-3),
                    ymax = mean_time + sd_time),
                width = 0.2, linewidth = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = method_colors) +
  labs(x = NULL, y = "Mean Runtime (seconds, log10)",
       title = "Computational Cost by Method") +
  theme_mda(legend_pos = "none")
ggsave("figs/timing_comparison.pdf", p, width = 20, height = 12, units = "cm")

# Coverage (MAR + MNAR, 5 parameters)
for (mech_i in c("MAR", "MNAR")) {
  for (pn in param_names) {
    fig <- cov_agg %>% filter(mech == mech_i, param == pn)
    p <- ggplot(fig, aes(x = miss, y = coverage, group = method)) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_line(aes(color = method), linewidth = 0.5) +
      geom_point(aes(color = method, shape = method), size = 1.0) +
      scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.8, 0.9, 0.95, 1)) +
      scale_method_aes +
      facet_grid(dist ~ N_label) +
      labs(x = "Missingness Rate", y = "95% CI Coverage",
           title = paste0("Coverage: ", param_labels[pn], " (", mech_i, ")")) +
      theme_mda()
    ggsave(sprintf("figs/coverage_%s_%s.pdf", tolower(mech_i), param_files[pn]), p,
           width = 30, height = 20, units = "cm")
  }
}

# MSE (MAR + MNAR, 5 parameters)
for (mech_i in c("MAR", "MNAR")) {
  for (pn in param_names) {
    fig <- mse_agg %>% filter(mech == mech_i, param == pn)
    p <- ggplot(fig, aes(x = miss, y = mse, group = method)) +
      geom_line(aes(color = method), linewidth = 0.5) +
      geom_point(aes(color = method, shape = method), size = 1.0) +
      scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
      scale_y_log10() +
      scale_method_aes +
      facet_grid(dist ~ N_label, scales = "free_y") +
      labs(x = "Missingness Rate", y = "MSE (log10)",
           title = paste0("MSE: ", param_labels[pn], " (", mech_i, ")")) +
      theme_mda()
    ggsave(sprintf("figs/mse_%s_%s.pdf", tolower(mech_i), param_files[pn]), p,
           width = 30, height = 20, units = "cm")
  }
}

cat(sprintf("\nDone. %d figures saved to figs/\n", total_figs))
