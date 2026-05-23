# helper: nearest positive-semidefinite projection
#
# Pairwise-deletion correlation/covariance matrices are not guaranteed to be
# positive semidefinite.  This function zeroes out negative eigenvalues and
# reconstructs, yielding the nearest PSD matrix (Higham 1988) in Frobenius
# norm.  A small ridge is *not* added; zeros are acceptable eigenvalues for a
# PSD matrix and are handled by the downstream validation in smriti_impute().
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
#' using a Lagrangian constraint.  Only originally-missing values are updated;
#' observed data is held fixed.
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param initial_imputation A matrix or data frame of the same dimensions as
#'   `data[, time_cols]`, containing initial imputed values. If `NULL`, a simple
#'   column-mean imputation is performed as a fallback.
#' @param lambda A numeric value specifying the per-observation penalty weight
#'   on the covariance-constraint term.  The covariance gradient is deliberately
#'   un-normalised (no division by n-1) so the constraint remains effective at
#'   any sample size.  Defaults to 1.0.  Increase for higher-dimensional or
#'   noisy targets; decrease for small samples.
#' @param learning_rate A numeric value for the gradient descent step size.
#'   Defaults to 0.001.
#' @param tol A numeric value for the internal convergence tolerance
#'   (Frobenius norm). Defaults to 1e-6.
#' @param max_iter An integer specifying the maximum number of iterations
#'   for the gradient descent projection. Defaults to 2000.
#' @param robust A logical value. Setting it to TRUE uses a robust target
#'   constructed from pairwise Spearman correlations and column-wise MAD,
#'   projected to the nearest positive-semidefinite matrix. This protects
#'   against outliers and heavy-tailed noise at the cost of some asymptotic
#'   efficiency under exact Gaussianity.
#' @param custom_target An optional p x p covariance matrix (where p is the
#'   number of longitudinal columns) to be used as the structural manifold.
#'   If provided, `robust` and `sigma_target` calculations are bypassed.
#'   This allows users to supply MAR-valid covariance targets from FIML or
#'   EM algorithms.
#'
#' @return A data frame with imputed and structurally refined values.
#'   Only the originally-missing cells are modified; observed values are
#'   returned unchanged.
#'
#' @examples
#' # Simulated longitudinal data with scattered missingness
#' df <- data.frame(
#'   T1 = c(1.2, NA, 2.8, 3.1),
#'   T2 = c(2.1, 2.5, NA, 4.0),
#'   T3 = c(3.0, 3.3, 4.1, NA)
#' )
#' smriti_impute(df, time_cols = 1:3, robust = FALSE)
#'
#' @export
smriti_impute <- function(data, time_cols, initial_imputation = NULL,
                          lambda = 1.0, learning_rate = 0.001, tol = 1e-6,
                          max_iter = 2000, robust = TRUE, custom_target = NULL) {

  # Type and dimension guards
  if (nrow(data) <= 1) {
    stop("Dataset must have >1 row to compute target covariance.")
  }
  if (!all(sapply(data[, time_cols], is.numeric))) {
    stop("All specified time_cols must be strictly numeric.")
  }

  # missingness mask
  x_raw   <- as.matrix(data[, time_cols])
  mask    <- ifelse(is.na(x_raw), 1.0, 0.0)
  storage.mode(mask) <- "double"   

  # check for entirely-missing columns
  na_counts <- colSums(is.na(data[, time_cols]))
  all_missing <- na_counts == nrow(data)
  if (any(all_missing)) {
    stop("Column(s) ", paste(names(which(all_missing)), collapse = ", "),
         " are 100% missing; target covariance cannot be estimated.")
  }

  # initial imputation
  if (is.null(initial_imputation)) {
    warning("initial_imputation is NULL. Falling back to simple column-mean ",
            "imputation. For better results, consider passing an initial ",
            "imputation from 'missRanger' or 'mice'.")
    
    x_hallucinated <- x_raw
    for (i in seq_len(ncol(x_hallucinated))) {
      na_idx <- is.na(x_hallucinated[, i])
      if (any(na_idx)) {
        x_hallucinated[na_idx, i] <- mean(x_hallucinated[, i], na.rm = TRUE)
      }
    }
  } else {
    x_hallucinated <- as.matrix(initial_imputation)
  }

  # Terminal guard against user-supplied NA matrices
  if (any(is.na(x_hallucinated))) {
    stop("The initial_imputation matrix still contains missing values (NAs). ",
         "The initial imputation must be complete before covariance projection.")
  }

  # target covariance from raw (incomplete) data
  if (!is.null(custom_target)) {
    sigma_target <- as.matrix(custom_target)
    if (nrow(sigma_target) != length(time_cols) || ncol(sigma_target) != length(time_cols)) {
      stop("Dimensions of custom_target must match the number of time_cols.")
    }
    sigma_target <- nearest_psd(sigma_target)
  } else if (robust) {
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

  # validate target conditioning
  # Allow zero eigenvalues (valid PSD matrix); reject only negative ones.
  target_eig <- eigen(sigma_target, symmetric = TRUE,
                      only.values = TRUE)$values
  if (any(target_eig < -1e-12)) {
    stop("Target covariance matrix is not positive semidefinite (smallest ",
         "eigenvalue = ", format(min(target_eig), digits = 3), "). ",
         "This should not happen after nearest_psd(). Please report as a bug.")
  }

  # C++ Lagrangian projection
  x_refined <- constrain_covariance(
    X_imp        = x_hallucinated,
    mask         = mask,
    Sigma_target = sigma_target,
    lambda       = lambda,
    lr           = learning_rate,
    max_iter     = max_iter,
    tol          = tol
  )

  # convergence diagnostic
  initial_cov  <- stats::cov(x_hallucinated)
  final_cov    <- stats::cov(x_refined)
  initial_dist <- sqrt(sum((initial_cov - sigma_target)^2))
  final_dist   <- sqrt(sum((final_cov   - sigma_target)^2))
  improvement  <- initial_dist - final_dist

  if (final_dist > tol) {
    warning(sprintf("Covariance projection did not reach tolerance. Final Dist: %.4e (Tol: %.4e).",
                    final_dist, tol))
  }

  # Guarantee original observed values remain perfectly untouched.
  # This protects against any scaling/drift in the initial_imputation matrix.
  obs_idx <- !is.na(x_raw)
  x_refined[obs_idx] <- x_raw[obs_idx]

  # verify observed data is untouched (post-injection guard)
  obs_after  <- x_refined[!is.na(x_raw)]
  if (max(abs(x_raw[!is.na(x_raw)] - obs_after)) > 1e-12) {
    warning("Observed values were unexpectedly modified during projection. ",
            "Maximum observed-value drift: ",
            format(max(abs(x_raw[!is.na(x_raw)] - obs_after)), digits = 3),
            ". This may indicate a bug in the masking logic.")
  }

  # assemble output
  final_data <- data
  final_data[, time_cols] <- x_refined
  final_data
}
