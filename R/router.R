#' Smriti Automated Longitudinal Imputation
#'
#' This function performs an automated routing and refinement for longitudinal
#' missing data. It establishes a target covariance manifold from observed data,
#' and then projects an initial imputation back toward the structural manifold 
#' using a Lagrangian constraint.
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param initial_imputation A matrix or data frame of the same dimensions as 
#'   `data[, time_cols]`, containing initial imputed values. If `NULL`, a simple 
#'   column-mean imputation is performed as a fallback.
#' @param lambda A numeric value specifying the penalty weight for the
#'   Lagrangian constraint. Defaults to 0.5.
#' @param learning_rate A numeric value for the gradient descent step size. 
#'   Defaults to 1e-7.
#' @param tol A numeric value for the convergence tolerance. Defaults to 1e-6.
#'   Note: This is currently reserved for future use in the C++ backend.
#' @param max_iter An integer specifying the maximum number of iterations. 
#'   Defaults to 1000.
#' @param robust A logical value. Setting it to TRUE uses a robust covariance 
#'   estimator (Spearman correlation and MAD) on the raw data to define the 
#'   target manifold, protecting against outliers and heavy-tailed skew.
#'
#' @return A data frame with imputed and structurally refined values.
#' @export
smriti_impute <- function(data, time_cols, initial_imputation = NULL, 
                          lambda = 0.5, learning_rate = 1e-7, tol = 1e-6, 
                          max_iter = 1000, robust = TRUE) {

  # Handle initial imputation fallback
  if (is.null(initial_imputation)) {
    warning("initial_imputation is NULL. Falling back to simple column-mean imputation. ",
            "For better results, consider passing an initial imputation from 'missRanger' or 'mice'.")
    
    x_hallucinated <- as.matrix(data[, time_cols])
    for (i in seq_len(ncol(x_hallucinated))) {
      na_idx <- is.na(x_hallucinated[, i])
      if (any(na_idx)) {
        x_hallucinated[na_idx, i] <- mean(x_hallucinated[, i], na.rm = TRUE)
      }
    }
  } else {
    x_hallucinated <- as.matrix(initial_imputation)
  }

  if (robust) {
    # Calculate robust target covariance from raw data (pairwise complete)
    # Using Spearman correlation and MAD as a lightweight robust proxy.
    # This avoids the catch-22 of calculating covariance on already-imputed data.
    cor_target <- stats::cor(data[, time_cols], use = "pairwise.complete.obs", method = "spearman")
    sd_target <- apply(data[, time_cols], 2, stats::mad, na.rm = TRUE)
    
    # Reconstruct covariance matrix from robust correlation and robust scales (D * R * D)
    sigma_target <- diag(sd_target) %*% cor_target %*% diag(sd_target)
  } else {
    # Establish the target manifold based on all available pairwise
    # information, prioritizing maximal use of observed data.
    sigma_target <- stats::cov(data[, time_cols], use = "pairwise.complete.obs")
  }

  # Call to C++ Armadillo backend
  x_refined <- constrain_covariance(
    X_imp = x_hallucinated,
    Sigma_target = sigma_target,
    lambda = lambda,
    lr = learning_rate,
    max_iter = max_iter
  )

  final_data <- data
  final_data[, time_cols] <- x_refined
  final_data
}
