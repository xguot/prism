library(ggplot2)
library(dplyr)
library(tidyr)

# load data if not already in environment
if (!exists("final_results")) {
        final_results <- readRDS("tests/simulation_results.rds")
}

plot_data <- final_results %>%
        pivot_longer(cols = c("FIML", "missForest", "Smriti_Nonrobust", "Smriti"), 
                     names_to = "Method", 
                     values_to = "RB")

# crop extreme outliers to prevent raw ml hallucinations from 
# compressing the scale of the smriti and fiml performance.
plot_data <- plot_data %>%
        filter(RB > -2, RB < 2)

# ensure logical ordering of methods for the x-axis
plot_data$Method <- factor(plot_data$Method, 
                           levels = c("FIML", "missForest", "Smriti_Nonrobust", "Smriti"))

ggplot(plot_data, aes(x = Method, y = RB, fill = Method)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        facet_grid(Dist ~ Miss, labeller = label_both) +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal() +
        labs(title = "Relative Bias of Slope Variance Estimate",
             subtitle = "Benchmarking Robustness-Efficiency Tradeoff",
             x = "Method",
             y = "Relative Bias (RB)") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("tests/simulation_results_plot.png", width = 12, height = 8)
