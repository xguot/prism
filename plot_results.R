suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Configuration
args     <- commandArgs(trailingOnly = TRUE)
rds_file <- if (length(args) > 0) args[1] else "sim_results/sim_prepped.rds"

stopifnot(file.exists(rds_file))

# Shared plot theme
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

# Shared scales and aesthetics
method_levels <- c("FIML", "TSRE", "BRITS", "CSDI", "GAIN", "SAITS", "NIMIWAE", "notMIWAE")

method_colors <- c(
  "FIML"       = "#E41A1C",
  "TSRE"       = "#FF7F00",
  "BRITS"      = "#377EB8",
  "CSDI"       = "#4DAF4A",
  "GAIN"       = "#984EA3",
  "SAITS"      = "#00BFC4",
  "NIMIWAE"    = "#F0E442",
  "notMIWAE"   = "#000000"
)

method_linetypes <- c(
  "FIML"       = "dashed",
  "TSRE"       = "dotdash",
  "BRITS"      = "solid",
  "CSDI"       = "solid",
  "GAIN"       = "solid",
  "SAITS"      = "solid",
  "NIMIWAE"    = "solid",
  "notMIWAE"   = "dotted"
)

scale_method_aes <- list(
  scale_color_manual(values = method_colors, breaks = method_levels),
  scale_linetype_manual(values = method_linetypes, breaks = method_levels),
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
)

dist_levels <- c("Normal", "Log-normal", "t", "Outlier")
miss_breaks <- c(0.05, 0.15, 0.30)
miss_labels <- c("5%", "15%", "30%")
N_levels    <- c("N = 100", "N = 200", "N = 500", "N = 1k", "N = 5k", "N = 10k")

# Load and preprocess data
cat(sprintf("Loading %s ...\n", rds_file))
res <- readRDS(rds_file)

res <- res %>% filter(method != "missForest")

res$method   <- factor(res$method, levels = method_levels)
res$dist     <- factor(res$dist, levels = dist_levels)
res$N_label  <- factor(res$N_label, levels = N_levels)
res$miss_pct <- factor(paste0(res$miss_rate * 100, "%"), levels = miss_labels)

# Aggregate across replicates
agg <- res %>%
  group_by(N, N_label, miss_rate, miss_pct, dist, mech, method, param) %>%
  summarise(
    across(c(mse, rel_bias, coverage),
           list(mean = ~ mean(.x, na.rm = TRUE),
                sd   = ~ sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

params <- sort(unique(agg$param))
param_short <- c(
  "Var Intercept" = "var_int",
  "Var Slope"     = "var_slope",
  "Covariance"    = "covar",
  "Intercept"     = "int",
  "Slope"         = "slope"
)

dir.create("figs", showWarnings = FALSE)

# Plot generation
for (pn in params) {
  ps <- param_short[pn]
  cat(sprintf("\nProcessing %s\n", pn))

  # Dynamic variables
  y_title <- ifelse(pn == "Covariance", "Absolute Raw Bias", "Absolute Relative Bias (%)")
  arb_threshold <- if (pn == "Covariance") {
    NULL
  } else {
    geom_hline(yintercept = 10, linetype = "dashed", color = "grey50", linewidth = 0.4)
  }

  # MSE MAR
  fig <- agg %>% filter(mech == "MAR", param == pn)
  p <- ggplot(fig, aes(x = miss_rate, y = mse_mean, group = method)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = "MSE", title = pn) +
    theme_mda()
  ggsave(sprintf("figs/mse_mar_%s.pdf", ps), p, width = 30, height = 20, units = "cm")

  # MSE MNAR
  fig <- agg %>% filter(mech == "MNAR", param == pn)
  p <- ggplot(fig, aes(x = miss_rate, y = mse_mean, group = method)) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = "MSE", title = pn) +
    theme_mda()
  ggsave(sprintf("figs/mse_mnar_%s.pdf", ps), p, width = 30, height = 20, units = "cm")

  # Absolute Relative / Raw Bias MAR
  fig <- agg %>% filter(mech == "MAR", param == pn)
  p <- ggplot(fig, aes(x = miss_rate, y = abs(rel_bias_mean), group = method)) +
    arb_threshold +  # Inject the dynamic threshold line
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = y_title, title = pn) + # Inject the dynamic label
    theme_mda()
  ggsave(sprintf("figs/arb_mar_%s.pdf", ps), p, width = 30, height = 20, units = "cm")

  # Absolute Relative / Raw Bias MNAR
  fig <- agg %>% filter(mech == "MNAR", param == pn)
  p <- ggplot(fig, aes(x = miss_rate, y = abs(rel_bias_mean), group = method)) +
    arb_threshold +  # Inject the dynamic threshold line
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_method_aes +
    facet_grid(dist ~ N_label, scales = "free_y") +
    labs(x = "Missingness Rate", y = y_title, title = pn) + # Inject the dynamic label
    theme_mda()
  ggsave(sprintf("figs/arb_mnar_%s.pdf", ps), p, width = 30, height = 20, units = "cm")

  # Coverage MAR
  fig <- agg %>% filter(mech == "MAR", param == pn)
  p <- ggplot(fig, aes(x = miss_rate, y = coverage_mean, group = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.50, 0.80, 0.90, 0.95, 1.0)) +
    scale_method_aes +
    facet_grid(dist ~ N_label) +
    labs(x = "Missingness Rate", y = "95% CI Coverage", title = pn) +
    theme_mda()
  ggsave(sprintf("figs/cov_mar_%s.pdf", ps), p, width = 30, height = 20, units = "cm")

  # Coverage MNAR
  fig <- agg %>% filter(mech == "MNAR", param == pn)
  p <- ggplot(fig, aes(x = miss_rate, y = coverage_mean, group = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(color = method), size = 1.0) +
    scale_x_continuous(breaks = miss_breaks, labels = miss_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.50, 0.80, 0.90, 0.95, 1.0)) +
    scale_method_aes +
    facet_grid(dist ~ N_label) +
    labs(x = "Missingness Rate", y = "95% CI Coverage", title = pn) +
    theme_mda()
  ggsave(sprintf("figs/cov_mnar_%s.pdf", ps), p, width = 30, height = 20, units = "cm")
}

cat(sprintf("\nExecution complete. %d figures saved to figs/\n", length(params) * 6))
