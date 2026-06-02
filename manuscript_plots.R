# ══════════════════════════════════════════════════════════════════════════════
# Manuscript-Level Performance Plots for smriti
# Uses existing prod_results.rds — no re-simulation needed.
# Plot conventions adapted from C1 simulation study (tutor reference).
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Shared plot theme (tutor convention) ─────────────────────────────────────
theme_smriti <- function() {
  theme_bw() +
    theme(
      legend.position      = c(0.85, 0.75),
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

# ══════════════════════════════════════════════════════════════════════════════
# Load & preprocess
# ══════════════════════════════════════════════════════════════════════════════
prod <- readRDS("sim_results/prod_results.rds")
prod$method   <- factor(prod$method, levels = method_levels)
prod$N_label  <- factor(paste0("N = ", prod$N),
                        levels = paste0("N = ", sort(unique(prod$N))))
prod$miss_pct <- factor(paste0(prod$miss * 100, "%"),
                        levels = c("5%", "15%", "30%"))

# ── Aggregate for summary plots ──────────────────────────────────────────────
agg <- prod %>%
  group_by(N, N_label, miss, miss_pct, dist, mech, method) %>%
  summarise(
    Frob     = mean(f_dist, na.rm = TRUE),
    Frob_sd  = sd(f_dist, na.rm = TRUE),
    Bias_pct = mean(s_var_bias, na.rm = TRUE),
    Bias_sd  = sd(s_var_bias, na.rm = TRUE),
    Time     = mean(time_sec, na.rm = TRUE),
    .groups  = "drop"
  )

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Frobenius Distance by Missingness Rate (MAR, headline plot)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 1: Frobenius Distance by Missingness Rate (MAR)...\n")

fig1_data <- agg %>% filter(mech == "MAR")

p1 <- ggplot(fig1_data, aes(x = miss, y = Frob, group = method)) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_continuous(
    breaks = c(0.05, 0.15, 0.30),
    labels = c("5%", "15%", "30%")
  ) +
  scale_method_aes +
  facet_grid(dist ~ N_label) +
  labs(
    x        = "Missingness Rate",
    y        = "Frobenius Distance to True Covariance",
    title    = "Covariance Recovery by Missingness Rate (MAR)",
    subtitle = "Lower is better.  Rows: distribution.  Columns: sample size."
  ) +
  theme_smriti()

ggsave("manuscript_figures/fig1_frobenius_by_miss.pdf",
       plot = p1, width = 27, height = 15, units = "cm")
ggsave("manuscript_figures/fig1_frobenius_by_miss.png",
       plot = p1, width = 27, height = 15, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Frobenius Distance by Sample Size (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 2: Frobenius Distance by Sample Size (MAR)...\n")

p2 <- ggplot(fig1_data, aes(x = N, y = Frob, group = method)) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000),
    labels = c("100", "200", "500", "1k", "5k")
  ) +
  scale_method_aes +
  facet_grid(dist ~ miss_pct) +
  labs(
    x        = "Sample Size (N)",
    y        = "Frobenius Distance to True Covariance",
    title    = "Covariance Recovery by Sample Size (MAR)",
    subtitle = "Lower is better.  Does advantage scale with N?"
  ) +
  theme_smriti()

ggsave("manuscript_figures/fig2_frobenius_by_N.pdf",
       plot = p2, width = 27, height = 15, units = "cm")
ggsave("manuscript_figures/fig2_frobenius_by_N.png",
       plot = p2, width = 27, height = 15, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Outlier Degradation (Δ Frobenius Normal → Outlier, MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 3: Outlier Degradation (MAR)...\n")

fig3_data <- agg %>%
  filter(mech == "MAR", dist %in% c("Normal", "Outlier")) %>%
  select(N, miss, dist, method, Frob) %>%
  pivot_wider(names_from = dist, values_from = Frob) %>%
  mutate(Delta = Outlier - Normal) %>%
  group_by(miss, method) %>%
  summarise(
    Delta_mean = mean(Delta, na.rm = TRUE),
    Delta_sd   = sd(Delta, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(miss_pct = factor(paste0(miss * 100, "%"), levels = c("5%", "15%", "30%")))

p3 <- ggplot(fig3_data, aes(x = method, y = Delta_mean, fill = method)) +
  geom_col(position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Delta_mean - Delta_sd, ymax = Delta_mean + Delta_sd),
                width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = method_colors, breaks = method_levels, guide = "none") +
  facet_wrap(~ miss_pct, nrow = 1) +
  labs(
    x        = "",
    y        = expression(Delta ~ "Frobenius (Outlier - Normal)"),
    title    = "Outlier-Induced Degradation in Covariance Recovery (MAR)",
    subtitle = expression(Lower ~ Delta * " = more robust to 5% contamination at each missingness rate.")
  ) +
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
# FIGURE 4 — Slope Variance Bias (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 4: Slope Variance Bias (MAR)...\n")

fig4_data <- agg %>% filter(mech == "MAR")

p4 <- ggplot(fig4_data, aes(x = N, y = Bias_pct, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.6) +
  geom_point(aes(color = method), size = 1.5) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000),
    labels = c("100", "200", "500", "1k", "5k")
  ) +
  scale_method_aes +
  facet_grid(dist ~ miss_pct) +
  labs(
    x        = "Sample Size (N)",
    y        = expression("Bias in Recovered " * sigma[s]^2 * " (%)"),
    title    = "Slope Variance Recovery by Sample Size (MAR)",
    subtitle = "Closer to zero is better.  Reference line at zero (unbiased)."
  ) +
  theme_smriti()

ggsave("manuscript_figures/fig4_slope_bias.pdf",
       plot = p4, width = 27, height = 15, units = "cm")
ggsave("manuscript_figures/fig4_slope_bias.png",
       plot = p4, width = 27, height = 15, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Timing Comparison (MAR, mean across conditions)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 5: Timing Comparison (MAR)...\n")

fig5_data <- agg %>%
  filter(mech == "MAR") %>%
  group_by(method) %>%
  summarise(
    Time_mean = mean(Time, na.rm = TRUE),
    Time_sd   = sd(Time, na.rm = TRUE),
    .groups   = "drop"
  )

p5 <- ggplot(fig5_data, aes(x = reorder(method, Time_mean), y = Time_mean, fill = method)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, Time_mean - Time_sd), ymax = Time_mean + Time_sd),
                width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = method_colors, breaks = method_levels, guide = "none") +
  labs(
    x        = "",
    y        = "Mean Time (seconds)",
    title    = "Computation Time by Method (MAR, pooled conditions)",
    subtitle = "Error bars: ± 1 SD.  Smriti routing step only (~0.6 s); full pipeline adds ~8 s."
  ) +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )

ggsave("manuscript_figures/fig5_timing.pdf",
       plot = p5, width = 18, height = 10, units = "cm")
ggsave("manuscript_figures/fig5_timing.png",
       plot = p5, width = 18, height = 10, units = "cm", dpi = 300)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Smriti_Robust vs Competitors: Frobenius at 30% Missingness (MAR)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 6: 30% Missingness Comparison (MAR)...\n")

fig6_data <- agg %>%
  filter(mech == "MAR", miss == 0.30)

p6 <- ggplot(fig6_data, aes(x = N, y = Frob, group = method)) +
  geom_line(aes(linetype = method, color = method), linewidth = 0.7) +
  geom_point(aes(color = method), size = 2.0) +
  scale_x_log10(
    breaks = c(100, 200, 500, 1000, 5000),
    labels = c("100", "200", "500", "1k", "5k")
  ) +
  scale_method_aes +
  facet_wrap(~ dist, nrow = 1) +
  labs(
    x        = "Sample Size (N)",
    y        = "Frobenius Distance to True Covariance",
    title    = "Covariance Recovery at 30% Missingness (MAR)",
    subtitle = "Does Smriti retain advantage under heavy missingness?"
  ) +
  theme_smriti()

ggsave("manuscript_figures/fig6_high_miss.pdf",
       plot = p6, width = 27, height = 10, units = "cm")
ggsave("manuscript_figures/fig6_high_miss.png",
       plot = p6, width = 27, height = 10, units = "cm", dpi = 300)

cat("\nAll 6 figures saved to manuscript_figures/\n")
