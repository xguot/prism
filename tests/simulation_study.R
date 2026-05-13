library(smriti)
library(missForest)
library(lavaan)

#' Generate Linear GCM Data with Latent Correlation
#' @param n Sample size.
#' @param dist Distribution type ("Normal" or "Lognormal").
#' @param t Number of time points (default 4).
generate_gcm_data <- function(n, dist = "Normal", t = 4) {
        rho <- 0.5
        Z_L <- rnorm(n)
        Z_S_raw <- rnorm(n)
        Z_S <- rho * Z_L + sqrt(1 - rho^2) * Z_S_raw

        if (dist == "Lognormal") {
                u_L <- (exp(Z_L) - exp(0.5)) / sqrt((exp(1) - 1) * exp(1))
                u_S <- (exp(Z_S) - exp(0.5)) / sqrt((exp(1) - 1) * exp(1))
        } else {
                u_L <- Z_L
                u_S <- Z_S
        }

        # Fixed effects: Intercept=6, Slope=2
        L <- 6 + u_L
        S <- 2 + u_S

        data <- matrix(0, nrow = n, ncol = t)
        for (i in 1:t) {
                # Standardize residual error
                e_raw <- rnorm(n)
                if (dist == "Lognormal") {
                        e <- (exp(e_raw) - exp(0.5)) / sqrt((exp(1) - 1) * exp(1))
                } else {
                        e <- e_raw
                }
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

#' Extract Slope Variance
#' @param fit A lavaan object.
get_slope_var <- function(fit) {
        if (!inherits(fit, "lavaan")) return(NA)
        pt <- parameterEstimates(fit)
        val <- pt$est[pt$lhs == "S" & pt$op == "~~" & pt$rhs == "S"]
        if (length(val) == 0) return(NA)
        return(val)
}

#' Run Simulation Experiment
#' @param reps Number of replications.
#' @param n Sample size.
#' @param miss_rate Missingness rate.
#' @param dist Distribution type.
run_experiment <- function(reps = 100, n = 200, miss_rate = 0.1, dist = "Normal") {
        results <- data.frame(FIML = numeric(reps), 
                              MissForest = numeric(reps),
                              Smriti_Nonrobust = numeric(reps),
                              Smriti = numeric(reps))
        
        gcm_model <- '
                L =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
                S =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
                L ~~ 0*S  # Misspecification: incorrectly fixing covariance to 0
                L ~~ L
                S ~~ S
        '

        for (i in 1:reps) {
                clean_data <- generate_gcm_data(n, dist = dist)
                miss_data <- introduce_missingness(clean_data, miss_rate)
                
                # 1. FIML
                fit_fiml <- try(growth(gcm_model, data = miss_data, missing = "fiml"), silent = TRUE)
                results$FIML[i] <- get_slope_var(fit_fiml)
                
                # 2. Raw missForest
                imp_mf <- try(missForest(miss_data)$ximp, silent = TRUE)
                if (is.data.frame(imp_mf)) {
                        fit_mf <- try(growth(gcm_model, data = imp_mf), silent = TRUE)
                        results$MissForest[i] <- get_slope_var(fit_mf)
                } else {
                        results$MissForest[i] <- NA
                }

                # 3. Smriti (Non-robust)
                # Acts as the high-efficiency baseline for Gaussian data.
                imp_sn <- try(smriti_impute(miss_data, time_cols = 1:4, robust = FALSE), silent = TRUE)
                if (is.data.frame(imp_sn)) {
                        fit_sn <- try(growth(gcm_model, data = imp_sn), silent = TRUE)
                        results$Smriti_Nonrobust[i] <- get_slope_var(fit_sn)
                } else {
                        results$Smriti_Nonrobust[i] <- NA
                }
                
                # 4. Smriti (Robust)
                imp_sm <- try(smriti_impute(miss_data, time_cols = 1:4), silent = TRUE)
                if (is.data.frame(imp_sm)) {
                        fit_sm <- try(growth(gcm_model, data = imp_sm), silent = TRUE)
                        results$Smriti[i] <- get_slope_var(fit_sm)
                } else {
                        results$Smriti[i] <- NA
                }
        }
        
        # True slope variance is 1.0
        rb_results <- as.data.frame(apply(results, 2, function(x) calc_rb(x, 1.0)))
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
