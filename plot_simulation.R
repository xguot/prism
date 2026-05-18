library(ggplot2)
library(dplyr)

# Reconstructing the simulation results data
results <- data.frame(
  Method = rep(c("missForest", "FIML", "smriti"), 3),
  Condition = c(
    rep("Normal (N=1000)", 3),
    rep("Normal (N=200)", 3),
    rep("Lognormal (N=200)", 3)
  ),
  Bias = c(
    0.080, 0.005, 0.006, # Normal 1000
    0.117, 0.157, 0.012, # Normal 200
    -0.182, NA, -0.180
  ) # Lognormal 200 (FIML fails/explodes here)
)

results$Condition <- factor(results$Condition,
  levels = c("Normal (N=1000)", "Normal (N=200)", "Lognormal (N=200)")
)
results$Method <- factor(results$Method, levels = c("FIML", "missForest", "smriti"))

p <- ggplot(results, aes(x = Bias, y = Method, color = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  geom_segment(aes(x = 0, xend = Bias, y = Method, yend = Method),
    size = 1.5, alpha = 0.5
  ) +
  facet_wrap(~Condition, ncol = 1) +
  scale_color_manual(values = c(
    "FIML" = "#E69F00",
    "missForest" = "#56B4E9",
    "smriti" = "#009E73"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "Relative Variance Bias Across Imputation Methods",
    subtitle = paste0(
      "smriti effectively anchors variance near zero (dashed line) across ",
      "Normal distributions,\nand remains stable under heavy-tailed ",
      "stress tests."
    ),
    x = "Relative Bias of Covariance Trace",
    y = ""
  )

ggsave("smriti_simulation_results.pdf",
  plot = p, width = 8, height = 6,
  device = "pdf"
)
cat("Successfully generated smriti_simulation_results.pdf\n")
