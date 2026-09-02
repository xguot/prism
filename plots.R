suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("load_results.R")
prod <- load_results()

# Output root: override for safe test runs (PRISM_FIG_DIR=/tmp/...)
fig_root <- Sys.getenv("PRISM_FIG_DIR", "figs")

# Shared theme — adapted from mda/plot_results.R (unchanged style)
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

# Method levels and aesthetics (unchanged)
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

# Factor levels (unchanged)
dist_levels <- c("Normal", "t5", "Outlier", "Lognormal")
miss_breaks <- c(0.05, 0.10, 0.15, 0.30)
miss_labels <- c("5%", "10%", "15%", "30%")

prod$method  <- factor(prod$method, levels = method_levels)
prod$dist    <- factor(prod$dist, levels = dist_levels)
prod$N_label <- factor(paste0("N = ", prod$N),
                       levels = paste0("N = ", c(100, 200, 500, 1000, 5000, 10000)))

# Methods kept in the main text
keep_methods <- c("FIML", "MICE (RF)", "missForest", "PRISM", "PRISM_MI")

# Per-model parameter definitions (GCM file names match the original figs/)
model_defs <- list(
  GCM = list(
    beta_true  = c(psi_L = 1, psi_S = 1, psi_LS = 0, beta_L = 6, beta_S = 2),
    labels     = c(psi_L = "Var Intercept", psi_S = "Var Slope",
                   psi_LS = "Covariance", beta_L = "Intercept", beta_S = "Slope"),
    files      = c(psi_L = "Var_Intercept", psi_S = "Var_Slope",
                   psi_LS = "Covariance", beta_L = "Intercept", beta_S = "Slope"),
    est_cols   = c(psi_L = "est_var_L", psi_S = "est_var_S",
                   psi_LS = "est_cov_LS", beta_L = "est_L", beta_S = "est_S"),
    se_cols    = c(psi_L = "se_var_L", psi_S = "se_var_S",
                   psi_LS = "se_cov_LS", beta_L = "se_L", beta_S = "se_S")
  ),
  SEM = list(
    beta_true  = c(b21 = 0.5, psi_F1 = 1, psi_F2 = 0.75),
    labels     = c(b21 = "Path eta2 ~ eta1", psi_F1 = "Var Factor 1",
                   psi_F2 = "Resid Var Factor 2"),
    files      = c(b21 = "Path_b21", psi_F1 = "Var_Factor1",
                   psi_F2 = "ResidVar_Factor2"),
    est_cols   = c(b21 = "est_b21", psi_F1 = "est_var_F1", psi_F2 = "est_var_F2"),
    se_cols    = c(b21 = "se_b21",  psi_F1 = "se_var_F1",  psi_F2 = "se_var_F2")
  )
)

total_figs <- 0L

for (model_i in c("GCM", "SEM")) {
  def <- model_defs[[model_i]]
  beta_true   <- def$beta_true
  param_names <- names(beta_true)
  param_labels <- def$labels
  param_files  <- def$files
  est_cols    <- def$est_cols
  se_cols     <- def$se_cols

  fig_dir <- file.path(fig_root, model_i)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("Plotting model %s -> %s/\n", model_i, fig_dir))
  d <- prod %>% filter(model == model_i)

  # Aggregated summaries use the median and MAD (SD-consistent), matching
  # the robust pooling in analysis.R: a few degenerate fits in the
  # heavy-tail cells otherwise dominate mean-based curves by orders of
  # magnitude.  Rates (coverage, convergence) keep their mean definition.

  # Aggregate Frobenius
  frob_agg <- d %>%
    filter(is.finite(f_dist), method %in% keep_methods) %>%
    group_by(N, N_label, miss, dist, mech, method) %>%
    summarise(frob_med = stats::median(f_dist, na.rm = TRUE),
              frob_mad = stats::mad(f_dist, na.rm = TRUE), .groups = "drop")

  # Per-parameter bias long table
  bias_agg <- do.call(rbind, lapply(param_names, function(pn) {
    tv <- beta_true[pn]; ec <- est_cols[pn]
    d %>%
      filter(is.finite(.data[[ec]]), method %in% keep_methods) %>%
      mutate(param = pn, est = .data[[ec]],
             rel_bias = if (abs(tv) < 1e-12) est - tv else 100 * (est - tv) / tv) %>%
      select(N, N_label, miss, dist, mech, method, param, rel_bias)
  })) %>%
    group_by(N, N_label, miss, dist, mech, method, param) %>%
    summarise(arb_med = stats::median(abs(rel_bias), na.rm = TRUE),
              arb_mad = stats::mad(abs(rel_bias), na.rm = TRUE), .groups = "drop")

  # Aggregate SE ratio
  se_ratio_agg <- do.call(rbind, lapply(param_names, function(pn) {
    sc <- se_cols[pn]; ec <- est_cols[pn]
    if (!sc %in% names(d)) return(NULL)
    d %>%
      filter(is.finite(.data[[sc]]), method %in% keep_methods) %>%
      group_by(N, N_label, miss, dist, mech, method) %>%
      summarise(
        empirical_mad  = stats::mad(.data[[ec]], na.rm = TRUE),
        median_model_se = stats::median(.data[[sc]], na.rm = TRUE),
        se_ratio       = median_model_se / empirical_mad,
        .groups = "drop"
      ) %>%
      mutate(param = pn)
  }))

  # Aggregate convergence rate
  conv_agg <- d %>%
    filter(method %in% keep_methods) %>%
    group_by(N, N_label, miss, dist, mech, method) %>%
    summarise(conv_rate = mean(converged, na.rm = TRUE), .groups = "drop")

  # Aggregate timing
  timing_agg <- d %>%
    filter(method %in% keep_methods, is.finite(time_sec)) %>%
    group_by(method) %>%
    summarise(median_time = stats::median(time_sec, na.rm = TRUE),
              mad_time    = stats::mad(time_sec, na.rm = TRUE), .groups = "drop")

  # Aggregate coverage (approx: est +/- 1.96 * se)
  cov_agg <- do.call(rbind, lapply(param_names, function(pn) {
    tv <- beta_true[pn]; ec <- est_cols[pn]; sc <- se_cols[pn]
    if (!sc %in% names(d)) return(NULL)
    d %>%
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

  # Aggregate MSE (bias^2 + empirical variance, robustly)
  mse_agg <- do.call(rbind, lapply(param_names, function(pn) {
    tv <- beta_true[pn]; ec <- est_cols[pn]
    d %>%
      filter(is.finite(.data[[ec]]), method %in% keep_methods) %>%
      group_by(N, N_label, miss, dist, mech, method) %>%
      summarise(
        mse = (stats::median(.data[[ec]], na.rm = TRUE) - tv)^2 +
              stats::mad(.data[[ec]], na.rm = TRUE)^2,
        .groups = "drop"
      ) %>%
      mutate(param = pn)
  }))

  # FIGURE 1 — Frobenius MAR
  p <- frob_agg %>% filter(mech == "MAR") %>%
    ggplot(aes(x = miss, y = frob_med, group = method)) +
    geom_line(aes(color = method), linewidth = 0.5) +
    geom_point(aes(color = method, shape = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_y_log10() +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = "Frobenius Distance (log10)", title = "Frobenius Distance to True Covariance (MAR)") +
    theme_mda()
  ggsave(file.path(fig_dir, "Frobenius_Distance_MAR.pdf"), p, width = 30, height = 20, units = "cm")

  # FIGURE 2 — Frobenius MNAR
  p <- frob_agg %>% filter(mech == "MNAR") %>%
    ggplot(aes(x = miss, y = frob_med, group = method)) +
    geom_line(aes(color = method), linewidth = 0.5) +
    geom_point(aes(color = method, shape = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_y_log10() +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = "Frobenius Distance (log10)", title = "Frobenius Distance to True Covariance (MNAR)") +
    theme_mda()
  ggsave(file.path(fig_dir, "Frobenius_Distance_MNAR.pdf"), p, width = 30, height = 20, units = "cm")
  total_figs <- total_figs + 2L

  # Absolute Relative Bias per parameter (MAR + MNAR)
  for (mech_i in c("MAR", "MNAR")) {
    for (pn in param_names) {
      cat(sprintf("Plotting Bias (%s): %s...\n", mech_i, pn))
      tv <- beta_true[pn]
      y_lab <- if (abs(tv) < 1e-12) "Absolute Raw Bias" else paste0("Absolute Relative Bias in ", param_labels[pn], " (%)")
      hline <- if (abs(tv) < 1e-12) NULL else
        geom_hline(yintercept = 10, linetype = "dashed", color = "grey50", linewidth = 0.4)

      fig <- bias_agg %>% filter(mech == mech_i, param == pn)
      p <- ggplot(fig, aes(x = miss, y = arb_med, group = method)) +
        hline +
        geom_line(aes(color = method), linewidth = 0.5) +
        geom_point(aes(color = method, shape = method), size = 1.0) +
        scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
        scale_method_aes +
        facet_grid(dist ~ N_label, scales = "free_y") +
        labs(x = "Missingness Rate", y = y_lab,
             title = paste0(param_labels[pn], " (", mech_i, ")")) +
        theme_mda()
      ggsave(file.path(fig_dir, sprintf("bias_%s_%s.pdf", tolower(mech_i), param_files[pn])), p,
             width = 30, height = 20, units = "cm")
      total_figs <- total_figs + 1L
    }
  }

  # SE metrics per parameter (MAR + MNAR)
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
          labs(x = "Missingness Rate", y = "SE Ratio (Median Model SE / Empirical MAD)",
               title = paste0("SE Ratio: ", param_labels[pn], " (", mech_i, ")")) +
          theme_mda()
        ggsave(file.path(fig_dir, sprintf("se_ratio_%s_%s.pdf", tolower(mech_i), param_files[pn])), p,
               width = 30, height = 20, units = "cm")

        # Empirical MAD
        p <- ggplot(fig, aes(x = miss, y = empirical_mad, group = method)) +
          geom_line(aes(color = method), linewidth = 0.5) +
          geom_point(aes(color = method, shape = method), size = 1.0) +
          scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
          scale_y_log10() +
          scale_method_aes +
          facet_grid(dist ~ N_label, scales = "free_y") +
          labs(x = "Missingness Rate", y = "Empirical MAD (log10)",
               title = paste0("Empirical MAD: ", param_labels[pn], " (", mech_i, ")")) +
          theme_mda()
        ggsave(file.path(fig_dir, sprintf("empirical_se_%s_%s.pdf", tolower(mech_i), param_files[pn])), p,
               width = 30, height = 20, units = "cm")

        # Median Model SE
        p <- ggplot(fig, aes(x = miss, y = median_model_se, group = method)) +
          geom_line(aes(color = method), linewidth = 0.5) +
          geom_point(aes(color = method, shape = method), size = 1.0) +
          scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
          scale_y_log10() +
          scale_method_aes +
          facet_grid(dist ~ N_label, scales = "free_y") +
          labs(x = "Missingness Rate", y = "Median Model SE (log10)",
               title = paste0("Median Model SE: ", param_labels[pn], " (", mech_i, ")")) +
          theme_mda()
        ggsave(file.path(fig_dir, sprintf("model_se_%s_%s.pdf", tolower(mech_i), param_files[pn])), p,
               width = 30, height = 20, units = "cm")
        total_figs <- total_figs + 3L
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
    ggsave(file.path(fig_dir, sprintf("conv_rate_%s.pdf", tolower(mech_i))), p,
           width = 30, height = 20, units = "cm")
    total_figs <- total_figs + 1L
  }

  # Timing comparison (bar)
  p <- timing_agg %>%
    ggplot(aes(x = method, y = median_time, fill = method)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(median_time - mad_time, 1e-3),
                      ymax = median_time + mad_time),
                  width = 0.2, linewidth = 0.3) +
    scale_y_log10() +
    scale_fill_manual(values = method_colors) +
    labs(x = NULL, y = "Median Runtime (seconds, log10)",
         title = "Computational Cost by Method") +
    theme_mda(legend_pos = "none")
  ggsave(file.path(fig_dir, "timing_comparison.pdf"), p, width = 20, height = 12, units = "cm")
  total_figs <- total_figs + 1L

  # Coverage (MAR + MNAR, all parameters)
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
      ggsave(file.path(fig_dir, sprintf("coverage_%s_%s.pdf", tolower(mech_i), param_files[pn])), p,
             width = 30, height = 20, units = "cm")
      total_figs <- total_figs + 1L
    }
  }

  # MSE (MAR + MNAR, all parameters)
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
      ggsave(file.path(fig_dir, sprintf("mse_%s_%s.pdf", tolower(mech_i), param_files[pn])), p,
             width = 30, height = 20, units = "cm")
      total_figs <- total_figs + 1L
    }
  }
}

cat(sprintf("\nDone. %d figures saved to %s/{GCM,SEM}/\n", total_figs, fig_root))
