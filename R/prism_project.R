# helper: nearest positive-semidefinite projection
#
# Pairwise-deletion correlation/covariance matrices are not guaranteed to be
# positive semidefinite.  This function zeroes out negative eigenvalues and
# reconstructs, yielding the nearest PSD matrix (Higham 1988) in Frobenius
# norm.  A small ridge is *not* added; zeros are acceptable eigenvalues for a
# PSD matrix and are handled by the downstream validation.
#' @keywords internal
nearest_psd <- function(mat) {
  if (any(is.na(mat))) {
    stop("Input matrix to nearest_psd contains NAs.")
  }
  eig <- eigen(mat, symmetric = TRUE)
  vals <- eig$values
  if (all(vals > 1e-12)) {
    return(mat)
  }
  vals[vals < 0] <- 0
  result <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  (result + t(result)) / 2
}


# Internal covariance projection engine.
#
# Project an initial imputation matrix onto the structural manifold defined
# by a supplied target covariance matrix, using Lagrangian-constrained
# gradient descent.  Only originally-missing cells are updated.
#
#' @keywords internal
prism_project <- function(data, time_cols, sigma_target,
                            initial_imputation = NULL,
                            lambda = 1.0, learning_rate = 0.001,
                            tol = 1e-6, max_iter = 2000) {

  # Type and dimension guards
  if (nrow(data) <= 1) {
    stop("Dataset must have >1 row to compute target covariance.")
  }
  if (!all(sapply(data[, time_cols], is.numeric))) {
    stop("All specified time_cols must be strictly numeric.")
  }

  # dimension and type validation for initial_imputation
  if (!is.null(initial_imputation)) {
    if (nrow(initial_imputation) != nrow(data) ||
        ncol(initial_imputation) != length(time_cols)) {
      stop("Dimensions of initial_imputation (",
           nrow(initial_imputation), "x", ncol(initial_imputation),
           ") must match the target data subset (",
           nrow(data), "x", length(time_cols), ").")
    }
    if (!all(sapply(initial_imputation, is.numeric))) {
      stop("The initial_imputation matrix must be strictly numeric.")
    }
  }

  # Validate sigma_target dimensions
  sigma_target <- as.matrix(sigma_target)
  if (nrow(sigma_target) != length(time_cols) ||
      ncol(sigma_target) != length(time_cols)) {
    stop("Dimensions of sigma_target must match the number of time_cols.")
  }
  sigma_target <- nearest_psd(sigma_target)

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

  # Terminal guard against user-supplied NA/Inf matrices
  if (any(is.na(x_hallucinated)) || any(is.infinite(x_hallucinated))) {
    stop("The initial_imputation matrix contains NAs or Infs. ",
         "The initial imputation must be complete and finite.")
  }

  # Pairwise completeness guard: check for NA in covariance target
  if (any(is.na(sigma_target))) {
    stop("Target covariance contains NAs.")
  }

  # validate target conditioning
  target_eig <- eigen(sigma_target, symmetric = TRUE,
                      only.values = TRUE)$values
  if (any(target_eig < -1e-12)) {
    stop("Target covariance matrix is not positive semidefinite (smallest ",
         "eigenvalue = ", format(min(target_eig), digits = 3), "). ",
         "This should not happen after nearest_psd(). Please report as a bug.")
  }

  # Standardize the data to prevent gradient blow-up in C++
  col_means <- apply(x_hallucinated, 2, mean)
  col_sds   <- apply(x_hallucinated, 2, stats::sd)

  col_sds[col_sds < 1e-10] <- 1.0

  x_scaled <- scale(x_hallucinated, center = col_means, scale = col_sds)

  # Scale the target covariance to match the scaled data
  scaling_mat <- diag(1/col_sds)
  sigma_scaled <- scaling_mat %*% sigma_target %*% scaling_mat

  # C++ Lagrangian projection
  cpp_result <- constrain_covariance(
    X_imp        = x_scaled,
    mask         = mask,
    Sigma_target = sigma_scaled,
    lambda       = lambda,
    lr           = learning_rate,
    max_iter     = max_iter,
    tol          = tol
  )
  x_refined_scaled <- cpp_result$X_refined

  # Final PSD enforcement on the refined covariance manifold
  sigma_refined_scaled <- stats::cov(x_refined_scaled)
  sigma_psd_scaled <- nearest_psd(sigma_refined_scaled)

  if (sqrt(sum((sigma_refined_scaled - sigma_psd_scaled)^2)) > 1e-9) {
    e_old <- eigen(sigma_refined_scaled, symmetric = TRUE)
    e_new <- eigen(sigma_psd_scaled, symmetric = TRUE)

    d_old_inv <- ifelse(e_old$values > 1e-12, 1/sqrt(e_old$values), 0)
    d_new     <- sqrt(e_new$values)

    s_old_inv_sqrt <- e_old$vectors %*% diag(d_old_inv) %*% t(e_old$vectors)
    s_new_sqrt     <- e_new$vectors %*% diag(d_new)     %*% t(e_new$vectors)

    x_centered <- x_refined_scaled -
      rep(colMeans(x_refined_scaled), each = nrow(x_refined_scaled))
    x_transformed <- (x_centered %*% s_old_inv_sqrt %*% s_new_sqrt) +
      rep(colMeans(x_refined_scaled), each = nrow(x_refined_scaled))

    x_refined_scaled <- x_refined_scaled +
      (x_transformed - x_refined_scaled) * mask
  }

  # Unscale the results back to original units
  x_refined <- t(t(x_refined_scaled) * col_sds + col_means)

  # convergence diagnostic
  initial_cov  <- stats::cov(x_hallucinated)
  final_cov    <- stats::cov(x_refined)
  initial_dist <- sqrt(sum((initial_cov - sigma_target)^2))
  final_dist   <- sqrt(sum((final_cov   - sigma_target)^2))

  if (!cpp_result$converged) {
    warning(sprintf(
      "C++ gradient descent did not converge in %d iterations. ",
      cpp_result$iterations
    ), "Final Frobenius distance to target: ",
    format(cpp_result$final_frob, digits = 3),
    ". Consider increasing max_iter or adjusting lambda.")
  }

  if (final_dist > tol) {
    warning(sprintf(
      "Covariance projection did not reach tolerance. ",
      "Final Dist: %.4e (Tol: %.4e).", final_dist, tol))
  }

  # Guarantee original observed values remain perfectly untouched
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
