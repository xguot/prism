library(smriti)
library(missForest)
library(lavaan)

#' Standardized Random Generation
#' @param n Sample size.
#' @param dist Distribution type ("Normal" or "Lognormal").
r_std <- function(n, dist = "Normal") {
        if (dist == "Normal") {
                return(rnorm(n, mean = 0, sd = 1))
        }

        /*
         * Shift and scale Lognormal(0,1) to Mean=0 and Var=1 to preserve
         * fixed effects while testing robustness to non-symmetric skew.
         */
        raw <- rlnorm(n, meanlog = 0, sdlog = 1)
        m <- exp(0.5)
        s <- sqrt((exp(1) - 1) * exp(1))
        return((raw - m) / s)
}

#' Generate Linear GCM Data
#' @param n Sample size.
#' @param dist Distribution type.
#' @param t Number of time points.
generate_gcm_data <- function(n, dist = "Normal", t = 4) {
        L_rand <- r_std(n, dist)
        S_rand <- r_std(n, dist)
        
        L <- 6 + L_rand
        S <- 2 + S_rand
        
        data <- matrix(0, nrow = n, ncol = t)
        for (i in 1:t) {
                e <- r_std(n, dist)
                data[, i] <- L + (i - 1) * S + e
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

#' Calculate Relative Bias
#' @param estimate Estimated value.
#' @param truth True value.
calc_rb <- function(estimate, truth) {
        return((estimate - truth) / truth)
}

#' Run Simulation Experiment
#' @param reps Number of replications.
#' @param n Sample size.
#' @param miss_rate Missingness rate.
#' @param dist Distribution type.
run_experiment <- function(reps = 100, n = 200, miss_rate = 0.1, dist = "Normal") {
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
                clean_data <- generate_gcm_data(n, dist = dist)
                miss_data <- introduce_missingness(clean_data, miss_rate)
                
                fit_fiml <- try(growth(gcm_model, data = miss_data, missing = "fiml"), silent = TRUE)
                results$FIML[i] <- if (inherits(fit_fiml, "lavaan")) coef(fit_fiml)["S~1"] else NA
                
                imp_mf <- try(missForest(miss_data)$ximp, silent = TRUE)
                if (is.data.frame(imp_mf)) {
                        fit_mf <- try(growth(gcm_model, data = imp_mf), silent = TRUE)
                        results$MissForest[i] <- if (inherits(fit_mf, "lavaan")) coef(fit_mf)["S~1"] else NA
                } else {
                        results$MissForest[i] <- NA
                }
                
                imp_sm <- try(smriti_impute(miss_data, time_cols = 1:4), silent = TRUE)
                if (is.data.frame(imp_sm)) {
                        fit_sm <- try(growth(gcm_model, data = imp_sm), silent = TRUE)
                        results$Smriti[i] <- if (inherits(fit_sm, "lavaan")) coef(fit_sm)["S~1"] else NA
                } else {
                        results$Smriti[i] <- NA
                }
        }
        
        rb_results <- as.data.frame(apply(results, 2, function(x) calc_rb(x, true_slope)))
        rb_results$N <- n
        rb_results$Miss <- miss_rate
        rb_results$Dist <- dist
        return(rb_results)
}

# Main Execution
conditions <- expand.grid(N = c(200, 1000), Miss = c(0.1, 0.3), Dist = c("Normal", "Lognormal"))
final_results_list <- list()

for (j in 1:nrow(conditions)) {
        cond <- conditions[j, ]
        cat(sprintf("Running N=%d, Miss=%.1f, Dist=%s\n", cond$N, cond$Miss, cond$Dist))
        final_results_list[[j]] <- run_experiment(reps = 100, n = cond$N, miss_rate = cond$Miss, dist = as.character(cond$Dist))
}

final_results <- do.call(rbind, final_results_list)
saveRDS(final_results, "tests/simulation_results.rds")
