#' Smriti Automated Longitudinal Imputation
#'
#' This function performs an automated routing and refinement for longitudinal
#' missing data. It establishes a target covariance manifold from observed data,
#' performs initial machine learning imputation, and then projects the result
#' back toward the structural manifold using a Lagrangian constraint.
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param lambda A numeric value specifying the penalty weight for the
#'   Lagrangian constraint.
#'
#' @importFrom stats shapiro.test cov
#' @importFrom missForest missForest
#' @export
smriti_impute <- function(data, time_cols, lambda = 0.1) {
        # 1. Diagnostic: Check if we should even use ML
        miss_rate <- sum(is.na(data)) / prod(dim(data))
        is_normal <- all(apply(data[, time_cols], 2,
                         function(x) stats::shapiro.test(x)$p.value > 0.05))

        # 2. Establish "Ground Truth" Covariance (The Target)
        # We use FIML or pairwise covariance on observed data
        sigma_target <- stats::cov(data[, time_cols], use = "pairwise.complete.obs")

        # 3. Raw Machine Learning Imputation
        # Using missForest to get the initial "hallucinated" values
        raw_imp_obj <- missForest::missForest(data)
        x_hallucinated <- as.matrix(raw_imp_obj$ximp[, time_cols])

        # 4. The Lagrangian Refinement (Your C++ Backend)
        # We project the ML output back toward the target manifold
        x_refined <- constrain_covariance(
                X_imp = x_hallucinated,
                Sigma_target = sigma_target,
                lambda = lambda,
                lr = 0.01,
                max_iter = 100
        )

        final_data <- data
        final_data[, time_cols] <- x_refined
        return(final_data)
}

#' Calculate Relative Bias
#'
#' @param estimates A numeric vector of parameter estimates.
#' @param true_val The true parameter value.
#'
#' @return The relative bias in percentage.
#' @export
calc_rb <- function(estimates, true_val) {
        100 * (mean(estimates, na.rm = TRUE) - true_val) / true_val
}
