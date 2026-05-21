library(smriti)

# Check 1: 100% Missing Column Edge Case
set.seed(20260521)
n <- 50
df <- data.frame(
  V1 = rnorm(n),
  V2 = rnorm(n),
  V3 = rep(NA_real_, n),
  V4 = rnorm(n)
)

cat("=== Check 1: 100% Missing Column Edge Case ===\n")
cat("Data dimensions:", nrow(df), "x", ncol(df), "\n")
cat("NA count V3:", sum(is.na(df$V3)), "/", nrow(df), "\n\n")

result <- tryCatch(
  smriti_impute(data = df, time_cols = c("V1", "V2", "V3", "V4"), robust = FALSE),
  error = function(e) {
    cat("ERROR CAUGHT:\n")
    cat("  Message:", conditionMessage(e), "\n")
    cat("  Class:", class(e)[1], "\n")
    
    msg <- conditionMessage(e)
    check1 <- grepl("100% missing", msg, ignore.case = FALSE)
    check2 <- grepl("V3", msg)
    check3 <- grepl("target covariance", msg)
    check4 <- !grepl("infinite", msg, ignore.case = TRUE)
    check5 <- !grepl("missing values.*eigen", msg, ignore.case = TRUE)
    check6 <- !grepl("missing.*infinite.*x", msg, ignore.case = TRUE)
    
    cat("\nVALIDATION:\n")
    cat("  [", if(check1) "PASS" else "FAIL", "] Message mentions '100% missing'\n")
    cat("  [", if(check2) "PASS" else "FAIL", "] Message identifies problematic column (V3)\n")
    cat("  [", if(check3) "PASS" else "FAIL", "] Message explains reason (target covariance)\n")
    cat("  [", if(check4) "PASS" else "FAIL", "] No 'infinite' error\n")
    cat("  [", if(check5) "PASS" else "FAIL", "] No eigen()-related crash\n")
    cat("  [", if(check6) "PASS" else "FAIL", "] No cryptic internal error\n")
    
    all_pass <- all(check1, check2, check3, check4, check5, check6)
    cat("\n  OVERALL: ", if(all_pass) "PASS" else "FAIL", "\n", sep = "")
    
    if (!all_pass) {
      cat("\n  Full message for diagnosis:\n  ", msg, "\n", sep = "")
    }
    
    quit(status = if(all_pass) 0 else 1)
  }
)

cat("FAIL: No error was raised for 100% missing column.\n")
quit(status = 1)
