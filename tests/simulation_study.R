library(smriti)
library(missForest)
library(lavaan)

#' Generate Linear GCM Data
#' @param n Sample size.
#' @param t Number of time points (default 4).
generate_gcm_data <- function(n, t = 4) {
        # Latent means: Intercept=6, Slope=2
        # Variances and Covariances = 1
        # Basis: 0, 1, 2, 3
        
        # Latent factors
        L <- rnorm(n, mean = 6, sd = 1)
        S <- rnorm(n, mean = 2, sd = 1)
        
        data <- matrix(0, nrow = n, ncol = t)
        for (i in 1:t) {
                # y = L + (i-1)*S + e
                data[, i] <- L + (i - 1) * S + rnorm(n, sd = 1)
        }
        
        colnames(data) <- paste0("T", 1:t)
        return(as.data.frame(data))
}

#' Introduce MCAR Missingness
#' @param data Data frame.
#' @param rate Missingness rate.
introduce_missingness <- function(data, rate) {
        n_cells <- prod(dim(data))
        indices <- sample(1:n_cells, size = floor(rate * n_cells))
        data_miss <- as.matrix(data)
        data_miss[indices] <- NA
        return(as.data.frame(data_miss))
}

#' Run Simulation Experiment
#' @param reps Number of replications.
#' @param n Sample size.
#' @param miss_rate Missingness rate.
run_experiment <- function(reps = 100, n = 200, miss_rate = 0.1) {
        true_slope <- 2
        results <- data.frame(FIML = numeric(reps), 
                              MissForest = numeric(reps), 
                              Smriti = numeric(reps))
        
        gcm_model <- '
                L =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
                S =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
                S ~ 1
        '

        for (i in 1:reps) {
                # Generate and break data
                clean_data <- generate_gcm_data(n)
                miss_data <- introduce_missingness(clean_data, miss_rate)
                
                # 1. FIML
                fit_fiml <- growth(gcm_model, data = miss_data, missing = "fiml")
                results$FIML[i] <- coef(fit_fiml)["S~1"]
                
                # 2. Raw missForest
                imp_mf <- missForest(miss_data)$ximp
                fit_mf <- growth(gcm_model, data = imp_mf)
                results$MissForest[i] <- coef(fit_mf)["S~1"]
                
                # 3. Smriti
                imp_smriti <- smriti_impute(miss_data, time_cols = 1:4)
                fit_smriti <- growth(gcm_model, data = imp_smriti)
                results$Smriti[i] <- coef(fit_smriti)["S~1"]
        }
        
        # Summarize Relative Bias
        summary <- colMeans(apply(results, 2, function(x) calc_rb(x, true_slope)))
        return(summary)
}

# Main Execution
conditions <- expand.grid(N = c(200, 1000), Miss = c(0.1, 0.3))
final_results <- list()

for (j in 1:nrow(conditions)) {
        cond <- conditions[j, ]
        cat(sprintf("Running N=%d, Miss=%.1f\n", cond$N, cond$Miss))
        final_results[[j]] <- run_experiment(reps = 100, n = cond$N, miss_rate = cond$Miss)
}

print(do.call(rbind, final_results))
