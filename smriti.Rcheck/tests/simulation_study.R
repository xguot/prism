library(smriti)
library(missForest)
library(lavaan)

# Command-line argument parsing for HPC scalability
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
        idx <- which(args == name)
        if (length(idx) > 0 && length(args) >= idx + 1) return(args[idx + 1])
        return(default)
}

REPS <- as.numeric(get_arg("--reps", if (identical(Sys.getenv("NOT_CRAN"), "true")) 100 else 1))
SEED <- as.numeric(get_arg("--seed", 42))

set.seed(SEED)

#' Construct a Latent Growth Curve Model (GCM) framework to simulate longitudinal trajectories.
#'
#' The underlying data-generating process incorporates latent correlation between
#' intercept and slope to represent realistic structural dependencies.
#'
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

        L <- 6 + u_L
        S <- 2 + u_S

        data <- matrix(0, nrow = n, ncol = t)
        for (i in 1:t) {
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


#' Introduce missingness following a Missing Completely At Random (MCAR) mechanism.
#'
#' @param data Data frame.
#' @param rate Missingness rate.
introduce_missingness <- function(data, rate) {
        n_cells <- prod(dim(data))
        indices <- sample(1:n_cells, size = floor(rate * n_cells))
        data_miss <- as.matrix(data)
        data_miss[indices] <- NA
        
        # Ensure no empty rows
        all_na <- apply(data_miss, 1, function(x) all(is.na(x)))
        if (any(all_na)) {
                for (i in which(all_na)) {
                        # Restore one random cell in the empty row
                        data_miss[i, sample(1:ncol(data_miss), 1)] <- data[i, sample(1:ncol(data), 1)]
                }
        }
        return(as.data.frame(data_miss))
}

#' Quantify the proportional deviation of the estimate from the population parameter.
#'
#' @param estimate Estimated value.
#' @param truth True value.
calc_rb <- function(estimate, truth) {
        return((estimate - truth) / truth)
}

#' Extract the slope variance parameter from the structural model estimates.
#'
#' @param fit A lavaan object.
get_slope_var <- function(fit) {
        if (!inherits(fit, "lavaan")) return(NA)
        pt <- parameterEstimates(fit)
        val <- pt$est[pt$lhs == "S" & pt$op == "~~" & pt$rhs == "S"]
        if (length(val) == 0) return(NA)
        return(val)
}

#' Execute a Monte Carlo simulation to evaluate recovery across diverse estimators.
#'
#' @param reps Number of replications.
#' @param n Sample size.
#' @param miss_rate Missingness rate.
#' @param dist Distribution type.
run_experiment <- function(reps = 100, n = 200, miss_rate = 0.1, dist = "Normal") {
        results <- data.frame(FIML = numeric(reps), 
                              missForest = numeric(reps),
                              smriti_nonrobust = numeric(reps),
                              smriti = numeric(reps))
        
        gcm_model <- '
                L =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
                S =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
                L ~~ 0*S  # Introduce structural misspecification to test robustness
                L ~~ L
                S ~~ S
        '

        for (i in 1:reps) {
                clean_data <- generate_gcm_data(n, dist = dist)
                miss_data <- introduce_missingness(clean_data, miss_rate)
                
                fit_fiml <- try(growth(gcm_model, data = miss_data, missing = "fiml"), silent = TRUE)
                results$FIML[i] <- get_slope_var(fit_fiml)
                
                imp_mf <- try(missForest(miss_data)$ximp, silent = TRUE)
                if (is.data.frame(imp_mf)) {
                        fit_mf <- try(growth(gcm_model, data = imp_mf), silent = TRUE)
                        results$missForest[i] <- get_slope_var(fit_mf)
                } else {
                        results$missForest[i] <- NA
                }

                # Evaluate the high-efficiency baseline for perfect Gaussian alignment
                imp_sn <- try(smriti_impute(miss_data, time_cols = 1:4, robust = FALSE), silent = TRUE)
                if (is.data.frame(imp_sn)) {
                        fit_sn <- try(growth(gcm_model, data = imp_sn), silent = TRUE)
                        results$smriti_nonrobust[i] <- get_slope_var(fit_sn)
                } else {
                        results$smriti_nonrobust[i] <- NA
                }
                
                imp_sm <- try(smriti_impute(miss_data, time_cols = 1:4), silent = TRUE)
                if (is.data.frame(imp_sm)) {
                        fit_sm <- try(growth(gcm_model, data = imp_sm), silent = TRUE)
                        results$smriti[i] <- get_slope_var(fit_sm)
                } else {
                        results$smriti[i] <- NA
                }
        }
        
        rb_results <- as.data.frame(apply(results, 2, function(x) calc_rb(x, 1.0)))
        rb_results$N <- n
        rb_results$Miss <- miss_rate
        rb_results$Dist <- dist
        return(rb_results)
}

# Orchestrate the simulation grid
if (sys.nframe() == 0) {
        conditions <- expand.grid(N = c(200, 1000), Miss = c(0.1, 0.3), Dist = c("Normal", "Lognormal"))
        final_results_list <- list()

        for (j in 1:nrow(conditions)) {
                cond <- conditions[j, ]
                cat(sprintf("Running N=%d, Miss=%.1f, Dist=%s with reps=%d\n", cond$N, cond$Miss, cond$Dist, REPS))
                final_results_list[[j]] <- run_experiment(reps = REPS, n = cond$N, miss_rate = cond$Miss, dist = as.character(cond$Dist))
        }

        final_results <- do.call(rbind, final_results_list)
        # Write to current directory; R CMD check runs tests in a temporary 'tests' folder
        saveRDS(final_results, "simulation_results.rds")
}
