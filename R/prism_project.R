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


# Emit a warning summarizing the optimizer termination state reported by the
# C++ engine.  Only "converged_feasible" is silent; every other state signals
# that the returned data must not be treated as a first-order solution to the
# unconstrained-within-feasible-set problem.
#' @keywords internal
warn_prism_diagnostics <- function(diag) {
  if (identical(diag$status, "converged_feasible")) {
    return(invisible(NULL))
  }
  msg <- switch(
    diag$status,
    constrained_geometric_limit = sprintf(
      paste0("Reached a constrained stationary point with covariance gap %.3e: ",
             "the mean-preservation constraints are binding. This is a genuine ",
             "geometric limit, not an early stop. Consider increasing ",
             "'lambda_sigma' or supplying a different initial imputation."),
      diag$feas_gap
    ),
    geometric_infeasible = sprintf(
      paste0("Reached a stationary point with inactive mean constraints but a ",
             "covariance gap of %.3e: the target is geometrically unreachable ",
             "given the frozen observed cells and the missingness pattern."),
      diag$feas_gap
    ),
    stalled_line_search = sprintf(
      paste0("Line search stalled before stationarity (KKT residual %.3e): ",
             "not a first-order solution. Check 'lr' and 'lambda_sigma'."),
      diag$r_kkT
    ),
    max_iter_reached = sprintf(
      paste0("Stopped at 'max_iter' before stationarity (KKT residual %.3e): ",
             "not a first-order solution. Increase 'max_iter' or adjust ",
             "'lambda_sigma'."),
      diag$r_kkT
    ),
    stop("Unknown optimizer status: ", diag$status, call. = FALSE)
  )
  warning(msg, call. = FALSE)
}


# Internal covariance projection engine (PRISM v2).
#
# Project an initial imputation matrix onto the structural target under the
# regularized dual objective
#
#   f(X) = 1/2 ||M o (X - X0)||_F^2 + (lambda_sigma / 2) ||S(X) - S_t||_F^2
#
# subject to the observed cells being frozen and the column means of the
# initial imputation preserved.  The C++ engine projects the gradient of the
# missing cells onto the zero-sum subspace of each column before every step
# (mean-preserving projected gradient descent) and reports KKT diagnostics:
# the projected-gradient norm (r_kkT), the covariance feasibility gap, and
# per-column Lagrange multiplier estimates.  Only originally-missing cells
# are updated.
#
# The data are standardized internally so that lambda_sigma and the
# tolerances are comparable across variable scales; results are unscaled
# afterwards and missing-cell column sums are re-anchored to the initial
# imputation as a defensive exactness guarantee.
#'
#' @keywords internal
prism_project <- function(data, ov_names, sigma_target,
                          initial_imputation = NULL,
                          lambda_sigma = NULL, lr = 0.01,
                          tol_kkT = 1e-4, tol_cov = 1e-6,
                          max_iter = 2000) {

  # structural weight defaults: the covariance term carries a 1/(n - 1)
  # factor, so the natural balance point scales with the sample size
  if (is.null(lambda_sigma)) {
    lambda_sigma <- nrow(data) / 2
  }

  # argument guards
  for (arg in c("lambda_sigma", "lr")) {
    val <- get(arg)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop("'", arg, "' must be a single finite number.", call. = FALSE)
    }
  }
  if (lambda_sigma < 0) {
    stop("'lambda_sigma' must be non-negative.", call. = FALSE)
  }
  if (lr <= 0) {
    stop("'lr' must be positive.", call. = FALSE)
  }
  for (arg in c("tol_kkT", "tol_cov")) {
    val <- get(arg)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0) {
      stop("'", arg, "' must be a single positive finite number.", call. = FALSE)
    }
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 ||
      max_iter < 1 || max_iter != round(max_iter)) {
    stop("'max_iter' must be a single positive integer.", call. = FALSE)
  }

  # Type and dimension guards
  if (nrow(data) <= 1) {
    stop("Dataset must have >1 row to compute target covariance.")
  }
  if (!all(vapply(data[, ov_names], is.numeric, logical(1)))) {
    stop("All observed variables must be strictly numeric.")
  }

  # dimension and type validation for initial_imputation
  if (!is.null(initial_imputation)) {
    if (nrow(initial_imputation) != nrow(data) ||
        ncol(initial_imputation) != length(ov_names)) {
      stop("Dimensions of initial_imputation (",
           nrow(initial_imputation), "x", ncol(initial_imputation),
           ") must match the observed variables of the data (",
           nrow(data), "x", length(ov_names), ").")
    }
    if (!all(vapply(initial_imputation, is.numeric, logical(1)))) {
      stop("The initial_imputation matrix must be strictly numeric.")
    }
  }

  # Validate sigma_target dimensions
  sigma_target <- as.matrix(sigma_target)
  if (nrow(sigma_target) != length(ov_names) ||
      ncol(sigma_target) != length(ov_names)) {
    stop("Dimensions of sigma_target must match the number of observed variables.")
  }
  sigma_target <- nearest_psd(sigma_target)

  # missingness mask
  x_raw <- as.matrix(data[, ov_names])
  mask  <- ifelse(is.na(x_raw), 1.0, 0.0)
  storage.mode(mask) <- "double"

  # check for entirely-missing columns
  na_counts <- colSums(is.na(data[, ov_names]))
  all_missing <- na_counts == nrow(data)
  if (any(all_missing)) {
    stop("Column(s) ", paste(names(which(all_missing)), collapse = ", "),
         " are 100% missing; target covariance cannot be estimated.")
  }

  # initial imputation
  if (is.null(initial_imputation)) {
    warning("initial_imputation is NULL. Falling back to simple column-mean ",
            "imputation. For better results, consider passing an initial ",
            "imputation from 'missRanger' or 'mice'.",
            call. = FALSE)

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

  # validate target conditioning
  target_eig <- eigen(sigma_target, symmetric = TRUE,
                      only.values = TRUE)$values
  if (any(target_eig < -1e-12)) {
    stop("Target covariance matrix is not positive semidefinite (smallest ",
         "eigenvalue = ", format(min(target_eig), digits = 3), "). ",
         "This should not happen after nearest_psd(). Please report as a bug.")
  }

  # Standardize so lambda_sigma and the tolerances are scale-invariant
  col_means <- apply(x_hallucinated, 2, mean)
  col_sds   <- apply(x_hallucinated, 2, stats::sd)

  col_sds[col_sds < 1e-10] <- 1.0

  x_scaled <- scale(x_hallucinated, center = col_means, scale = col_sds)

  # Scale the target covariance to match the scaled data
  scaling_mat  <- diag(1 / col_sds)
  sigma_scaled <- scaling_mat %*% sigma_target %*% scaling_mat

  # C++ v2 engine: regularized projection with mean-preserving gradients
  cpp_result <- constrain_covariance_v2(
    X_imp        = x_scaled,
    mask         = mask,
    Sigma_target = sigma_scaled,
    lambda_sigma = lambda_sigma,
    lr           = lr,
    max_iter     = max_iter,
    tol_kkT      = tol_kkT,
    tol_cov      = tol_cov
  )
  x_refined_scaled <- cpp_result$X_refined

  # Final PSD enforcement on the refined covariance manifold (defensive; the
  # sample covariance is a Gram matrix and PSD by construction)
  sigma_refined_scaled <- stats::cov(x_refined_scaled)
  sigma_psd_scaled <- nearest_psd(sigma_refined_scaled)

  if (sqrt(sum((sigma_refined_scaled - sigma_psd_scaled)^2)) > 1e-9) {
    e_old <- eigen(sigma_refined_scaled, symmetric = TRUE)
    e_new <- eigen(sigma_psd_scaled, symmetric = TRUE)

    d_old_inv <- ifelse(e_old$values > 1e-12, 1 / sqrt(e_old$values), 0)
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

  # Re-anchor missing-cell column sums to the initial imputation so that the
  # column means of X0 are preserved exactly in the returned data, even after
  # the defensive PSD fix-up and rescaling round-trips
  n_mis <- colSums(mask)
  if (any(n_mis > 0)) {
    anchor_sums   <- colSums(x_hallucinated * mask)
    current_sums  <- colSums(x_refined * mask)
    anchor_cols   <- which(n_mis > 0)
    for (j in anchor_cols) {
      x_refined[mask[, j] == 1, j] <- x_refined[mask[, j] == 1, j] +
        (anchor_sums[j] - current_sums[j]) / n_mis[j]
    }
  }

  # Guarantee original observed values remain perfectly untouched
  obs_idx <- !is.na(x_raw)
  x_refined[obs_idx] <- x_raw[obs_idx]

  # verify observed data is untouched (post-injection guard)
  obs_after <- x_refined[!is.na(x_raw)]
  if (max(abs(x_raw[!is.na(x_raw)] - obs_after)) > 1e-12) {
    warning("Observed values were unexpectedly modified during projection. ",
            "Maximum observed-value drift: ",
            format(max(abs(x_raw[!is.na(x_raw)] - obs_after)), digits = 3),
            ". This may indicate a bug in the masking logic.",
            call. = FALSE)
  }

  # optimization diagnostics attached to the completed data
  diag <- list(
    status       = cpp_result$status,
    converged    = cpp_result$converged,
    iterations   = cpp_result$iterations,
    r_kkT        = cpp_result$r_kkT,
    feas_gap     = cpp_result$feas_gap,
    objective    = cpp_result$objective,
    fidelity     = cpp_result$fidelity,
    cov_term     = cpp_result$cov_term,
    nu           = cpp_result$nu,
    grad_spread  = cpp_result$grad_spread,
    max_abs_nu   = cpp_result$max_abs_nu,
    lambda_sigma = lambda_sigma
  )
  warn_prism_diagnostics(diag)

  # assemble output
  final_data <- data
  final_data[, ov_names] <- x_refined
  attr(final_data, "prism_diagnostics") <- diag
  final_data
}
