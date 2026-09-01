#' Proper Multiple Imputation via Two-Level Stochastic Draws
#'
#' @description Generate \code{m} completed datasets by combining two levels
#'   of uncertainty, then projecting each pair onto its perturbed structural
#'   manifold:
#'   \enumerate{
#'     \item \strong{Missing-value (predictive) uncertainty}: each completed
#'           dataset starts from its own stochastic initial imputation
#'           \eqn{X^{(0,k)}}. By default a bootstrap-weighted missRanger
#'           forest is refit per draw; any stochastic imputer can be supplied
#'           through \code{initializer}.
#'     \item \strong{Parameter uncertainty}: free FIML parameter estimates
#'           \eqn{\hat{\theta}} and their asymptotic covariance matrix
#'           \eqn{\text{ACOV}(\hat{\theta})} are used to draw perturbed
#'           parameter vectors \eqn{\theta^{(k)} \sim
#'           \mathcal{N}(\hat{\theta}, \text{ACOV}(\hat{\theta}))}, from
#'           which perturbed model-implied covariances \eqn{\Sigma_k} are
#'           reconstructed.
#'   }
#'   Both sources of randomness enter the between-imputation variance under
#'   Rubin's rules, producing proper multiple imputations (Rubin 1987).
#'
#' @param data A data frame containing missing values.
#' @param fit A fitted lavaan object (must be estimated with
#'   \code{missing = "fiml"}).
#' @param m An integer giving the number of imputations. Must be at least 2.
#'   Defaults to 20.
#' @param initial_imputation Deprecated. A single deterministic initial
#'   imputation, reused for every draw; this disables level-1 uncertainty.
#'   Supply an \code{initializer} function instead.
#' @param initializer A function accepting \code{data} and returning a
#'   complete data frame or matrix for the observed variables, called once
#'   per imputation draw. It must be stochastic (e.g. a bootstrap- or
#'   predictive-draw-based imputer). If \code{NULL} (default), a
#'   bootstrap-weighted \code{missRanger} forest is used; the number of trees
#'   is controlled by \code{options("prism.num_trees")} (default 500).
#' @param seed A single integer seed for reproducibility of the complete draw
#'   sequence. If \code{NULL} (default), draws use the current RNG state.
#' @param return_mids A logical value. If \code{TRUE}, the result is cast to
#'   a \code{\link[mice]{mids}} object for compatibility with the
#'   \code{mice} pooling workflow. Requires the \code{mice} package.
#'   Defaults to \code{FALSE}.
#' @param target_means A logical or numeric vector controlling the mean
#'   target of each draw's projection. `TRUE` (default) targets the
#'   perturbed model-implied means of the observed variables (joint mean +
#'   covariance projection; the level-2 parameter draw perturbs the means
#'   as well). `FALSE` preserves the column means of each draw's own
#'   initial imputation exactly (legacy behaviour). A numeric vector of
#'   length p supplies a custom mean target. Columns without missing values
#'   keep their observed means in all modes.
#' @param ... Additional arguments passed to the internal projection engine,
#'   e.g. \code{lambda_sigma}, \code{lr}, \code{tol_kkT}, \code{tol_cov},
#'   \code{max_iter}.
#'
#' @details
#' The default initializer injects missing-value uncertainty through a
#' bootstrap-weighted forest: each call draws fresh Bayesian-bootstrap case
#' weights (\eqn{w_i \sim \text{Exp}(1)}) and refits \code{missRanger} on the
#' weighted sample. This varies the fitted forest across draws; it does not
#' add node-level predictive draws. For cell-level predictive uncertainty,
#' pass a custom \code{initializer} (e.g. one wrapping \code{miceRanger}
#' predictive matching). Within each draw, the PRISM v2 projection targets
#' the perturbed model-implied means by default (\code{target_means = TRUE}),
#' so the level-2 parameter draws propagate into the first moments as well;
#' \code{target_means = FALSE} preserves the column means of that draw's own
#' \eqn{X^{(0,k)}} instead. Parameter draws are sampled with eigenvalue
#' flooring and validated per draw; a draw whose implied moments cannot be
#' instantiated (degenerate fits in heavy-tail conditions) is re-rolled in a
#' shrinking neighbourhood of the estimates and, as a last resort, replaced
#' by the observed fit's implied moments with a warning.
#'
#' @return A list of \code{m} completed data frames (class
#'   \code{"prism_mi_list"}), each carrying a \code{prism_diagnostics}
#'   attribute, or a \code{mids} object if \code{return_mids = TRUE}.
#' @export
#'
#' @examples
#' \dontrun{
#' library(lavaan)
#' model <- "i =~ 1*T1 + 1*T2 + 1*T3
#'           s =~ 0*T1 + 1*T2 + 2*T3"
#' df <- data.frame(
#'   T1 = c(1.2, NA, 2.8, 3.1),
#'   T2 = c(2.1, 2.5, NA, 4.0),
#'   T3 = c(3.0, 3.3, 4.1, NA)
#' )
#' fit <- growth(model, data = df, missing = "fiml")
#' mi_list <- prism_mi(df, fit, m = 20, seed = 42)
#' }
prism_mi <- function(data, fit, m = 20, initial_imputation = NULL,
                     initializer = NULL, seed = NULL,
                     return_mids = FALSE, target_means = TRUE, ...) {
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.", call. = FALSE)
  }
  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      m < 2 || m != round(m)) {
    stop("'m' must be a single integer of at least 2.", call. = FALSE)
  }
  m <- as.integer(m)

  if (!is.null(initial_imputation) && !is.null(initializer)) {
    stop("Supply at most one of 'initial_imputation' or 'initializer'.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || seed != round(seed)) {
      stop("'seed' must be a single integer.", call. = FALSE)
    }
    set.seed(seed)
  }

  ov_names <- lavaan::lavNames(fit, "ov")

  # level 1: resolve the stochastic initializer
  deterministic_level1 <- FALSE
  if (is.null(initializer)) {
    if (!is.null(initial_imputation)) {
      warning(
        paste0("'initial_imputation' is deprecated and provides a ",
               "deterministic level-1 draw; missing-value uncertainty is ",
               "not propagated. Supply an 'initializer' function for ",
               "proper two-level multiple imputation."),
        call. = FALSE
      )
      x0_deterministic <- initial_imputation
      initializer <- function(data) x0_deterministic
      deterministic_level1 <- TRUE
    } else {
      num_trees <- getOption("prism.num_trees", 500L)
      initializer <- stochastic_rf_initializer(ov_names, num_trees)
    }
  }
  if (!is.function(initializer)) {
    stop("'initializer' must be a function accepting the data and returning ",
         "a complete data frame or matrix.", call. = FALSE)
  }

  theta_hat <- tryCatch(lavaan::coef(fit), error = function(e) NULL)
  acov <- tryCatch(as.matrix(lavaan::vcov(fit)), error = function(e) NULL)

  if (is.null(theta_hat) || is.null(acov) ||
      anyNA(theta_hat) || anyNA(acov)) {
    if (deterministic_level1) {
      stop("The fitted model's parameter estimates or asymptotic covariance ",
           "are unavailable (the model may not have converged) and the ",
           "initializer is deterministic, so no source of between-imputation ",
           "uncertainty is available. Refit the model or supply a ",
           "stochastic 'initializer'.", call. = FALSE)
    }
    warning("The fitted model's parameter estimates or asymptotic covariance ",
            "are unavailable (the model may not have converged); ",
            "parameter-uncertainty draws are disabled and each imputation ",
            "projects onto the observed fit's implied moments. ",
            "Between-imputation variance now comes from the stochastic ",
            "initializer only.", call. = FALSE)
    acov <- NULL
    draws <- NULL
  } else {
    # level 2: parameter draws from the asymptotic FIML sampling
    # distribution, with eigenvalue flooring so degenerate fits cannot
    # produce wild vectors
    draws <- prism_draw_parameters(theta_hat, acov, m)
  }

  pt <- lavaan::parTable(fit)
  free_idx <- which(pt$free > 0L)

  sample_cov  <- lavaan::lavInspect(fit, "sampstat")$cov
  sample_nobs <- lavaan::lavInspect(fit, "nobs")

  imputations <- vector("list", m)
  draw_fallbacks <- 0L
  extreme_fallbacks <- 0L

  for (i in seq_len(m)) {
    # level 1: fresh stochastic initial imputation
    x0_i <- validate_initializer_output(initializer(data), data, ov_names)

    # level 2: perturbed model-implied covariance from theta^(k)
    #
    # Even floored draws can yield non-finite or indefinite implied moments
    # in heavy-tail conditions; re-roll the draw in a shrinking
    # neighbourhood of the estimates (variance acov / attempt) and fall
    # back to the observed fit's implied moments as a last resort.
    fit_i <- NULL
    sigma_i <- NULL
    if (is.null(acov)) {
      # parameter draws are disabled: project every draw onto the observed
      # fit's implied moments when computable
      base_moments <- tryCatch(
        list(sigma = as.matrix(lavaan::fitted(fit)$cov), fit_i = fit),
        error = function(e) NULL
      )
      if (is.null(base_moments)) {
        # extreme fallback: the initial imputation's sample covariance
        # (between-imputation variance then comes from level 1 only)
        sigma_i <- stats::cov(x0_i)
        fit_i <- NULL
        extreme_fallbacks <- extreme_fallbacks + 1L
      } else {
        sigma_i <- base_moments$sigma
        fit_i <- base_moments$fit_i
        draw_fallbacks <- draw_fallbacks + 1L
      }
    } else {
      for (attempt in seq_len(10L)) {
        pt_i <- pt
        if (attempt == 1L) {
          pt_i$start[free_idx] <- draws[i, ]
          pt_i$est[free_idx]   <- draws[i, ]
        } else {
          d <- as.numeric(prism_draw_parameters(theta_hat, acov / attempt, 1L))
          pt_i$start[free_idx] <- d
          pt_i$est[free_idx]   <- d
        }

        # Suppress lavaan startup banner per iteration.  The model is not
        # refit (do.fit = FALSE): it is only instantiated to compute the
        # model-implied moments from the perturbed parameters, so the
        # post-fit convergence check must be disabled.
        fit_i <- tryCatch(
          suppressMessages(
            lavaan::lavaan(
              model      = pt_i,
              sample.cov = sample_cov,
              sample.nobs = sample_nobs,
              do.fit     = FALSE,
              se         = "none",
              test       = "none",
              check.post = FALSE
            )
          ),
          error = function(e) NULL
        )
        if (!is.null(fit_i)) {
          sigma_i <- tryCatch(as.matrix(lavaan::fitted(fit_i)$cov),
                              error = function(e) NULL)
          if (!is.null(sigma_i) && all(is.finite(sigma_i)) &&
              min(eigen(sigma_i, symmetric = TRUE,
                        only.values = TRUE)$values) > -1e-8) {
            break
          }
        }
        fit_i <- NULL
        sigma_i <- NULL
      }
      if (is.null(sigma_i)) {
        # fall back to the observed fit's implied moments when they can be
        # computed; a non-converged base fit can make even fitted() fail
        base_moments <- tryCatch(
          list(sigma = as.matrix(lavaan::fitted(fit)$cov), fit_i = fit),
          error = function(e) NULL
        )
        if (is.null(base_moments)) {
          sigma_i <- stats::cov(x0_i)
          fit_i <- NULL
          extreme_fallbacks <- extreme_fallbacks + 1L
        } else {
          sigma_i <- base_moments$sigma
          fit_i <- base_moments$fit_i
          draw_fallbacks <- draw_fallbacks + 1L
        }
      }
    }
    colnames(sigma_i) <- ov_names
    rownames(sigma_i) <- ov_names

    # Resolve the mean target for this draw from the perturbed fit (NULL
    # when the extreme fallback was used: resolve_mean_target then anchors
    # the initial imputation's means with a warning)
    mu_i <- tryCatch({
      mu <- lavaan::fitted(fit_i)$mean
      if (is.null(mu) || length(mu) != length(ov_names)) NULL
      else as.numeric(mu)
    }, error = function(e) NULL)
    mu_target <- resolve_mean_target(target_means, mu_i, ov_names)

    imputations[[i]] <- prism_project(
      data               = data,
      ov_names           = ov_names,
      sigma_target       = sigma_i,
      mu_target          = mu_target,
      initial_imputation = x0_i,
      ...
    )
    attr(imputations[[i]], "imputation") <- i
  }

  if (draw_fallbacks > 0L || extreme_fallbacks > 0L) {
    warning(
      sprintf("Parameter draws were degenerate in %d of %d imputations ",
              draw_fallbacks + extreme_fallbacks, m),
      sprintf("(%d used the observed fit's implied moments, %d used the ",
              draw_fallbacks, extreme_fallbacks),
      "initial imputation's sample covariance).",
      call. = FALSE
    )
  }

  class(imputations) <- c("prism_mi_list", "list")

  if (return_mids) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' is required for return_mids = TRUE. ",
           "Install it with install.packages('mice').", call. = FALSE)
    }
    # as.mids expects long format: the original incomplete data at
    # .imp = 0 followed by each imputation with its own .imp index
    orig <- data
    orig$.imp <- 0L
    long <- vector("list", m + 1L)
    long[[1L]] <- orig
    for (i in seq_len(m)) {
      d <- imputations[[i]]
      d$.imp <- i
      long[[i + 1L]] <- d
    }
    long <- do.call(rbind, long)
    return(suppressMessages(mice::as.mids(long)))
  }

  imputations
}


# Default level-1 initializer: a bootstrap-weighted missForest draw.
#
# Each call draws fresh Bayesian-bootstrap case weights and refits missRanger
# on the weighted sample, so successive calls return different imputations
# and propagate fitted-model uncertainty.  Node-level predictive draws are
# not performed; use a custom initializer (e.g. miceRanger) for that.
#' @keywords internal
stochastic_rf_initializer <- function(ov_names, num_trees) {
  force(ov_names)
  force(num_trees)
  function(data) {
    if (!requireNamespace("missRanger", quietly = TRUE)) {
      stop(
        paste0("The default stochastic initializer requires the ",
               "'missRanger' package. Install it or supply a custom ",
               "'initializer' function."),
        call. = FALSE
      )
    }
    sub <- data[, ov_names, drop = FALSE]
    w   <- stats::rexp(nrow(sub))
    missRanger::missRanger(
      sub,
      pmm.k       = 0,
      num.trees   = num_trees,
      case.weights = w,
      verbose     = 0
    )
  }
}


# Validate an initializer draw: complete, finite, correctly shaped, and not
# silently overwriting observed values.
#' @keywords internal
validate_initializer_output <- function(x0, data, ov_names) {
  if (is.data.frame(x0) &&
      all(ov_names %in% colnames(x0))) {
    x0 <- x0[, ov_names, drop = FALSE]
  } else if (!is.data.frame(x0) && !is.matrix(x0)) {
    stop("'initializer' must return a data frame or matrix.", call. = FALSE)
  }
  x0_mat <- as.matrix(x0)
  if (nrow(x0_mat) != nrow(data) ||
      ncol(x0_mat) != length(ov_names)) {
    stop("'initializer' must return an object with one row per observation ",
         "and one column per observed variable.", call. = FALSE)
  }
  if (!is.numeric(x0_mat)) {
    stop("'initializer' must return strictly numeric values.", call. = FALSE)
  }
  if (anyNA(x0_mat) || any(!is.finite(x0_mat))) {
    stop("'initializer' returned NA or non-finite values.", call. = FALSE)
  }

  obs_mask <- !is.na(as.matrix(data[, ov_names]))
  if (any(obs_mask) &&
      max(abs(x0_mat[obs_mask] -
              as.matrix(data[, ov_names])[obs_mask])) > 1e-6) {
    warning(
      paste0("The initializer modified observed values; they will be ",
             "reset to the original data during projection."),
      call. = FALSE
    )
  }
  x0_mat
}


# Draw perturbed parameter vectors theta ~ N(theta_hat, acov) with
# eigenvalue flooring.
#
# The asymptotic covariance of a degenerate fit can be near-singular or
# indefinite, and sampling from it directly produces wild parameter
# vectors whose implied moments are unusable (non-finite or non-PSD).
# Negative or near-zero eigenvalues are floored at a small multiple of the
# largest eigenvalue before drawing, which keeps the draws in a regular
# neighbourhood of the estimates without changing the procedure for
# well-conditioned fits.  The draw is exact for positive semidefinite acov
# whose smallest eigenvalue exceeds the floor.
#' @keywords internal
prism_draw_parameters <- function(theta_hat, acov, m) {
  eig <- eigen(acov, symmetric = TRUE)
  lam <- eig$values
  scale <- max(abs(lam), 1e-8)
  floor_val <- max(1e-10 * scale, 1e-8)
  lam_stable <- pmax(lam, floor_val)
  z <- matrix(stats::rnorm(m * length(theta_hat)),
              nrow = m, ncol = length(theta_hat))
  half <- eig$vectors %*% diag(sqrt(lam_stable)) %*% t(eig$vectors)
  matrix(theta_hat, nrow = m, ncol = length(theta_hat), byrow = TRUE) +
    z %*% half
}
