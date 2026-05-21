library(smriti)

# Case 1: All columns missing in one row
df_row_na <- data.frame(T1=c(1, NA, 3), T2=c(2, NA, 4))
cat("Testing all-missing row...\n")
tryCatch({
  res <- smriti_impute(df_row_na, time_cols=1:2)
  print(res)
}, error = function(e) cat("Error in all-missing row:", e$message, "\n"))

# Case 2: One column 100% missing
df_col_na <- data.frame(T1=c(1, 2, 3), T2=c(NA, NA, NA))
cat("\nTesting 100% missing column...\n")
tryCatch({
  res <- smriti_impute(df_col_na, time_cols=1:2)
  print(res)
}, error = function(e) cat("Error in 100% missing column:", e$message, "\n"))
