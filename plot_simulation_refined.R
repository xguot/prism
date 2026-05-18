library(ggplot2)
library(dplyr)

# Data reconstruction
results <- data.frame(
  Method = rep(c("FIML", "missForest", "smriti"), 3),
  Condition = c(
    rep("Normal (N=1000)", 3),
    rep("Normal (N=200)", 3),
    rep("Lognormal (N=200)", 3)
  ),
  Bias = c(
    0.005, 0.080, 0.006,  # Normal 1000
    0.157, 0.117, 0.012,  # Normal 200
    NA, -0.182, -0.180
  ) # Lognormal 200 (FIML NA)
)

results$Condition <- factor(results$Condition,
  levels = c("Normal (N=1000)", "Normal (N=200)", "Lognormal (N=200)")
)
results$Method <- factor(results$Method, levels = c("smriti", "missForest", "FIML"))

# Define colors
method_colors <- c(
  "smriti" = "#009E73",
  "missForest" = "#56B4E9",
  "FIML" = "#E69F00"
)

p <- ggplot(results, aes(x = Bias, y = Method, color = Method)) +
  # Reference line at zero
  geom_vline(
    xintercept = 0, linetype = "dashed", color = "gray50",
    linewidth = 0.8
  ) +
  # Lollipop lines
  geom_segment(aes(x = 0, xend = Bias, y = Method, yend = Method),
    linewidth = 1.2, alpha = 0.6
  ) +
  # Points
  geom_point(size = 4, alpha = 0.9) +
  # Faceting
  facet_wrap(~Condition, ncol = 1, strip.position = "left") +
  # Scales
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    # Align title and subtitle to the left baseline
    plot.title.position = "plot",
    plot.caption.position = "plot",

    # Text alignments and styles
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(
      size = 11, color = "gray30", hjust = 0,
      margin = margin(b = 15)
    ),

    # Legend at the top, horizontal
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.title = element_blank(),

    # Remove individual facet labels/names (we rely on legend)
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),

    # Grid lines
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),

    # Facet labels (Condition) on the left
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10, hjust = 1),
    strip.placement = "outside",

    # Spacing
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(
    title = "Imputation Method Variance Bias",
    subtitle = paste0(
      "SMRITI Anchoring: Consistent minimal variance bias in Normal tests, ",
      "highly stable with Lognormal."
    ),
    x = "Relative Bias of Covariance Trace",
    y = ""
  )

ggsave("smriti_simulation_refined.pdf",
  plot = p, width = 8, height = 7,
  device = "pdf"
)
ggsave("smriti_simulation_refined.png", plot = p, width = 8, height = 7, dpi = 300)
cat(paste0(
  "Successfully generated smriti_simulation_refined.pdf and ",
  "smriti_simulation_refined.png\n"
))
