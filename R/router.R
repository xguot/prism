# ---- helper: nearest positive-semidefinite projection ----
#
# Pairwise-deletion correlation/covariance matrices are not guaranteed to be
# positive semidefinite.  This function zeroes out negative eigenvalues and
# reconstructs, yielding the nearest PSD matrix (Higham 1988) in Frobenius
# norm.  A small ridge is *not* added — the resulting matrix may have zero
# eigenvalues, which is acceptable for the subsequent eigendecomposition
# validation but should trigger a user-facing error in smriti_impute().
#' @keywords internal
nearest_psd <- function(mat) {
  eig <- eigen(mat, symmetric = TRUE)
  vals <- eig$values
  if (all(vals > 1e-12)) {
    return(mat)                     # already PSD
  }
  vals[vals < 0] <- 0
  result <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  (result + t(result)) / 2          # symmetrise away floating-point drift
}


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
#'   covariance constraint relative to fidelity to the initial imputation.
#'   Higher values enforce the target covariance more strictly but allow
#'   greater deviation from the initial imputed values. Defaults to 1.0.
#' @param learning_rate A numeric value for the gradient descent step size.
#'   Defaults to 0.001.
#' @param tol A numeric value for warning threshold on final covariance
#'   deviation (Frobenius norm). Defaults to 1e-6.
#' @param max_iter An integer specifying the maximum number of iterations
#'   for the gradient descent projection. Defaults to 2000.
#' @param robust A logical value. Setting it to TRUE uses a robust target
#'   constructed from pairwise Spearman correlations and column-wise MAD,
#'   projected to the nearest positive-semidefinite matrix. This protects
#'   against outliers and heavy-tailed noise at the cost of some asymptotic
#'   efficiency under exact Gaussianity.
#'
#' @return A data frame with imputed and structurally refined values.
#' @export
smriti_impute <- function(data, time_cols, initial_imputation = NULL,
                          lambda = 1.0, learning_rate = 0.001, tol = 1e-6,
                          max_iter = 2000, robust = TRUE) {

  # ---- 1.  initial imputation ----
  if (is.null(initial_imputation)) {
    warning("initial_imputation is NULL. Falling back to simple column-mean ",
            "imputation. For better results, consider passing an initial ",
            "imputation from 'missRanger' or 'mice'.")

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

  # ---- 2.  target covariance from raw (incomplete) data ----
  if (robust) {
    # Robust path: pairwise Spearman correlation -> nearest PSD -> scale by
    # column MAD.  The PSD projection fixes the non-positive-definiteness
    # that pairwise deletion commonly introduces.
    cor_target <- stats::cor(data[, time_cols],
                             use = "pairwise.complete.obs",
                             method = "spearman")
    cor_target <- nearest_psd(cor_target)

    sd_target <- apply(data[, time_cols], 2, stats::mad, na.rm = TRUE)

    # D * R * D — now with a guaranteed-positive-semidefinite R
    sigma_target <- diag(sd_target) %*% cor_target %*% diag(sd_target)
    sigma_target <- nearest_psd(sigma_target)
  } else {
    # Non-robust path: pairwise Pearson covariance.  PSD projection is still
    # applied because pairwise deletion can produce indefinite matrices.
    sigma_target <- stats::cov(data[, time_cols],
                               use = "pairwise.complete.obs")
    sigma_target <- nearest_psd(sigma_target)
  }

  # ---- 3.  validate target conditioning ----
  target_eig <- eigen(sigma_target, symmetric = TRUE,
                      only.values = TRUE)$values
  if (min(target_eig) < 1e-12) {
    stop("Target covariance matrix is numerically singular (smallest ",
         "eigenvalue = ", format(min(target_eig), digits = 3), "). ",
         "The data may be too sparse across the specified time columns ",
         "to establish a valid covariance manifold. Consider reducing ",
         "the number of time columns or increasing the sample size.")
  }

  # ---- 4.  C++ Lagrangian projection ----
  x_refined <- constrain_covariance(
    X_imp        = x_hallucinated,
    Sigma_target = sigma_target,
    lambda       = lambda,
    lr           = learning_rate,
    max_iter     = max_iter
  )

  # ---- 5.  convergence diagnostic ----
  initial_cov  <- stats::cov(x_hallucinated)
  final_cov    <- stats::cov(x_refined)
  initial_dist <- sqrt(sum((initial_cov - sigma_target)^2))
  final_dist   <- sqrt(sum((final_cov   - sigma_target)^2))
  improvement  <- initial_dist - final_dist

  # Warn only if the projection made the fit *worse* (divergence), not
  # merely because it fell short of a tight tolerance on an approximate target.
  if (improvement < 0 && final_dist > tol) {
    warning(sprintf(
      paste0("Covariance projection diverged from target. ",
             "Final Frobenius distance: %.4e (initial: %.4e). ",
             "Consider reducing learning_rate (currently %.2e) ",
             "or verifying the target matrix conditioning."),
      final_dist, initial_dist, learning_rate
    ))
  }

  # ---- 6.  assemble output ----
  final_data <- data
  final_data[, time_cols] <- x_refined
  final_data
}
