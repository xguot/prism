# helper: nearest positive-semidefinite projection
#
# Pairwise-deletion correlation/covariance matrices are not guaranteed to be
# positive semidefinite.  This function zeroes out negative eigenvalues and
# reconstructs, yielding the nearest PSD matrix (Higham 1988) in Frobenius
# norm.  A small ridge is *not* added; zeros are acceptable eigenvalues for a
# PSD matrix and are handled by the downstream validation.  A materially
# non-PSD input (Frobenius change above numerical noise) triggers a warning
# rather than a silent retargeting.
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
  result <- (result + t(result)) / 2
  delta <- sqrt(sum((mat - result)^2))
  if (delta > 1e-6) {
    warning(
      "The covariance matrix was projected onto the PSD cone with Frobenius ",
      "change ", format(delta, digits = 3),
      ": this exceeds numerical noise and may indicate a misspecified ",
      "target (e.g. a non-converged model fit).",
      call. = FALSE
    )
  }
  result
}


# Emit a warning summarizing the optimizer termination state reported by the
# C++ engine.  Only "converged_feasible" is silent; "converged" with a
# nonzero covariance gap reflects the fidelity/structure trade-off of the
# regularized objective (not target infeasibility), and non-stationary stops
# are flagged explicitly.
#' @keywords internal
warn_prism_diagnostics <- function(diag) {
  if (identical(diag$status, "converged_feasible")) {
    return(invisible(NULL))
  }
  msg <- switch(
    diag$status,
    converged = sprintf(
      paste0("Reached a first-order stationary point with covariance gap ",
             "%.3e (tol_cov = %.3e): the regularized optimum does not fully ",
             "match the target. Increase 'lambda_sigma' to enforce the ",
             "model-implied structure more strongly."),
      diag$feas_gap, diag$tol_cov
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


# Resolve the mean target from the target_means argument.
#
# TRUE  -> model-implied means supplied by the caller (the default: joint
#          mean + covariance projection, consistent under MAR)
# FALSE -> NULL (the engine anchors the initial imputation's column means)
# numeric vector of length p -> custom mean target
# anything else (wrong length, non-numeric, non-logical, non-finite) -> error
#' @keywords internal
resolve_mean_target <- function(target_means, implied_means, ov_names) {
  p <- length(ov_names)
  if (is.numeric(target_means)) {
    if (length(target_means) != p) {
      stop("Length of 'target_means' (", length(target_means),
           ") must match the number of observed variables (", p, ").",
           call. = FALSE)
    }
    mu <- as.numeric(target_means)
    if (any(!is.finite(mu))) {
      stop("'target_means' must contain only finite values.", call. = FALSE)
    }
    return(mu)
  }
  if (is.logical(target_means) && length(target_means) == 1L &&
      !is.na(target_means)) {
    if (!target_means) {
      return(NULL)
    }
    if (is.null(implied_means) || length(implied_means) != p) {
      warning("Could not extract model-implied means from the lavaan fit; ",
              "falling back to mean preservation of the initial imputation.",
              call. = FALSE)
      return(NULL)
    }
    return(as.numeric(implied_means))
  }
  stop("'target_means' must be TRUE, FALSE, or a numeric vector of length ",
       p, ".", call. = FALSE)
}


# Internal covariance projection engine (PRISM v2).
#
# Project an initial imputation matrix onto the structural target under the
# regularized constrained objective
#
#   f(X) = 1 / (2 N_mis) ||M o (X - X0)||_F^2
#        + lambda_sigma / (2 p^2) ||S(X) - S_t||_F^2
#
# where N_mis is the number of missing cells and p the number of observed
# variables.  Both terms are normalized (fidelity per missing cell,
# covariance discrepancy per matrix entry), so lambda_sigma is dimensionless
# and comparable across sample sizes, variable counts, and missingness
# rates.  The objective is minimized subject to the observed cells being
# frozen and each column mean being fixed: to the scaled mean target when
# mu_target is supplied (a joint (mu, Sigma) projection), or to the initial
# imputation's column mean otherwise (the legacy anchor).  The C++ engine
# projects the gradient of the missing cells onto the zero-sum subspace of
# each column before every step and reports KKT diagnostics: the
# projected-gradient norm (r_kkT), the covariance feasibility gap, and
# per-column Lagrange multiplier estimates.  Only originally-missing cells
# are updated, and the returned matrix is the certified constrained optimum
# the diagnostics describe.
#
# The data are standardized internally so that the objective and tolerances
# are comparable across variable scales; results are unscaled afterwards and
# (in the legacy X0-anchor mode only) missing-cell column sums are
# re-anchored to the initial imputation as a defensive exactness guarantee.
#'
#' @keywords internal
prism_project <- function(data, ov_names, sigma_target,
                          mu_target = NULL,
                          initial_imputation = NULL,
                          lambda_sigma = NULL, lr = 1,
                          tol_kkT = 1e-4, tol_cov = 1e-6,
                          max_iter = 2000) {

  # structural weight default: with normalized losses lambda_sigma is
  # dimensionless; 1 gives equal weight to both terms (a sensitivity grid of
  # 0.25-4 is recommended for simulation studies)
  if (is.null(lambda_sigma)) {
    lambda_sigma <- 1.0
  }

  # mu_target validation
  if (!is.null(mu_target)) {
    mu_target <- as.numeric(mu_target)
    if (length(mu_target) != length(ov_names)) {
      stop("Length of 'mu_target' (", length(mu_target),
           ") must match the number of observed variables (",
           length(ov_names), ").", call. = FALSE)
    }
    if (any(!is.finite(mu_target))) {
      stop("'mu_target' must contain only finite values.", call. = FALSE)
    }
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

  if (sum(mask) == 0) {
    warning("No missing values in the observed variables; nothing to impute.",
            call. = FALSE)
  }

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

  # Enforce the observed cells to the original data before standardization
  # and anchor computation, so the preserved column means refer to the
  # completed data as it is actually constructed
  obs_cells <- !is.na(x_raw)
  if (any(obs_cells)) {
    dev <- max(abs(x_hallucinated[obs_cells] - x_raw[obs_cells]))
    if (dev > 1e-6) {
      warning("The initial imputation modified observed values (max ",
              "deviation ", format(dev, digits = 3),
              "); they are reset to the original data.", call. = FALSE)
    }
    x_hallucinated[obs_cells] <- x_raw[obs_cells]
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

  # Scale the mean target into the engine's standardized space:
  # mu_scaled_j = (mu_target_j - mean(X0[,j])) / sd(X0[,j]).  The engine
  # anchors the missing-cell column sums to n * mu_scaled_j - obs_sum_j, so
  # the completed column mean equals mu_target_j exactly (the affine
  # unscaling preserves this in the original units).  NULL keeps the legacy
  # X0-anchor behaviour.
  mu_scaled <- NULL
  if (!is.null(mu_target)) {
    mu_scaled <- (mu_target - col_means) / col_sds
  }

  # C++ v2 engine: regularized projection with mean anchoring; the returned
  # matrix is the certified constrained optimum described by the diagnostics
  cpp_result <- constrain_covariance_v2(
    X_imp        = x_scaled,
    mask         = mask,
    Sigma_target = sigma_scaled,
    lambda_sigma = lambda_sigma,
    lr           = lr,
    max_iter     = max_iter,
    tol_kkT      = tol_kkT,
    tol_cov      = tol_cov,
    mu_scaled_target = mu_scaled
  )
  x_refined_scaled <- cpp_result$X_refined

  # Defensive guard: the sample covariance is a Gram matrix and positive
  # semidefinite by construction, so a materially negative eigenvalue would
  # indicate a computational problem rather than something to repair
  eig_min <- min(eigen(stats::cov(x_refined_scaled), symmetric = TRUE,
                       only.values = TRUE)$values)
  if (eig_min < -1e-8) {
    warning("The refined covariance has a negative eigenvalue (",
            format(eig_min, digits = 3),
            ") below the numerical-noise floor; this suggests a ",
            "computational problem in the projection engine.",
            call. = FALSE)
  }

  # Unscale the results back to original units
  x_refined <- t(t(x_refined_scaled) * col_sds + col_means)

  # ── Mean handling ──────────────────────────────────────────────────────────
  # The engine already enforces the mean anchor internally (scaled mean
  # target when mu_target is supplied, X0 anchor otherwise), so the returned
  # matrix is the certified constrained optimum and no post-hoc correction
  # is applied in mean-targeting mode.  The legacy re-anchor below is a
  # defensive no-op for the X0-anchor mode only; it guards against drift in
  # the rescaling round-trips.
  n_mis <- colSums(mask)
  mean_gap <- NA_real_

  if (!is.null(mu_target) && any(n_mis > 0)) {
    # Mean-to-target gap of the RETURNED matrix, restricted to columns with
    # missing cells (the only columns the engine can move).  Fully observed
    # columns keep their observed means and cannot reach a target that
    # differs from them.
    mean_gap <- max(abs(colMeans(x_refined)[n_mis > 0] -
                        mu_target[n_mis > 0]))
  } else if (any(n_mis > 0)) {
    # Re-anchor missing-cell column sums to the initial imputation so that
    # the column means of X0 are preserved exactly (defensive no-op; the
    # engine already preserves them)
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
    n_mis        = cpp_result$n_mis,
    lambda_sigma = lambda_sigma,
    tol_cov      = tol_cov,
    mean_gap     = mean_gap,
    mu_target    = mu_target
  )
  warn_prism_diagnostics(diag)

  # assemble output
  final_data <- data
  final_data[, ov_names] <- x_refined
  attr(final_data, "prism_diagnostics") <- diag
  final_data
}
