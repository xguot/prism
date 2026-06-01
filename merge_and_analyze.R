library(dplyr)

# List all individual condition results
result_files <- list.files("sim_results", pattern = "prod_results_\\d+\\.rds", full.names = TRUE)

cat("Found", length(result_files), "individual result files.\n")

if (length(result_files) > 0) {
  # Read and combine all files
  all_data <- lapply(result_files, readRDS) %>% bind_rows()
  
  # Save consolidated file
  saveRDS(all_data, "sim_results/prod_results.rds")
  cat("Successfully merged all results into sim_results/prod_results.rds\n")
  
  # Perform analysis
  summary_stats <- all_data %>%
    group_by(method, dist, mech) %>%
    summarize(
      mean_bias = mean(s_var_bias, na.rm = TRUE),
      sd_s_var = sd(s_var, na.rm = TRUE),
      mean_time = mean(time_sec, na.rm = TRUE),
      n_reps = n(),
      .groups = 'drop'
    )
    
  print(summary_stats, n = 100)
} else {
  cat("No result files found in sim_results/ matching the pattern.\n")
}
