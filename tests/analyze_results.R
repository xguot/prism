library(ggplot2)
library(dplyr)
library(tidyr)

# Load data if not already in environment
if (!exists("final_results")) {
        final_results <- readRDS("tests/simulation_results.rds")
}

plot_data <- final_results %>%
        pivot_longer(cols = c(FIML, MissForest, Smriti), 
                     names_to = "Method", 
                     values_to = "RB")

# Crop extreme outliers to prevent raw ML hallucinations from 
# compressing the scale of the Smriti and FIML performance.
plot_data <- plot_data %>%
        filter(RB > -2, RB < 2)

ggplot(plot_data, aes(x = Method, y = RB, fill = Method)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        facet_grid(Dist ~ Miss, labeller = label_both) +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal() +
        labs(title = "Relative Bias of Latent Slope Estimate",
             subtitle = "Comparison of FIML, missForest, and Smriti across conditions",
             x = "Method",
             y = "Relative Bias (RB)") +
        theme(legend.position = "none")

ggsave("tests/simulation_results_plot.png", width = 10, height = 8)
