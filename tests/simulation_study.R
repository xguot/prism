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

#' Calculate Relative Bias
#' @param estimate Estimated value.
#' @param truth True value.
calc_rb <- function(estimate, truth) {
        return((estimate - truth) / truth)
}
