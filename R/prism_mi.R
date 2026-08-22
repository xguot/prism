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
#' predictive matching). Within each draw, the PRISM v2 projection preserves
#' the column means of that draw's own \eqn{X^{(0,k)}}, so mean uncertainty
#' propagates into the between-imputation variance as well.
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
                     return_mids = FALSE, ...) {
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
    } else {
      num_trees <- getOption("prism.num_trees", 500L)
      initializer <- stochastic_rf_initializer(ov_names, num_trees)
    }
  }
  if (!is.function(initializer)) {
    stop("'initializer' must be a function accepting the data and returning ",
         "a complete data frame or matrix.", call. = FALSE)
  }

  # level 2: parameter draws from the asymptotic FIML sampling distribution
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for prism_mi(). ",
         "Install it with install.packages('MASS').", call. = FALSE)
  }

  theta_hat <- lavaan::coef(fit)
  acov <- as.matrix(lavaan::vcov(fit))

  if (anyNA(theta_hat) || anyNA(acov)) {
    stop("Parameter estimates or asymptotic covariance matrix contain NAs. ",
         "The model may not be identified.", call. = FALSE)
  }

  draws <- tryCatch(
    MASS::mvrnorm(n = m, mu = theta_hat, Sigma = acov),
    error = function(e) {
      # Fall back to perturbing each parameter independently
      matrix(stats::rnorm(m * length(theta_hat), mean = theta_hat,
                          sd = sqrt(pmax(diag(acov), 1e-8))),
             nrow = m, byrow = TRUE)
    }
  )

  pt <- lavaan::parTable(fit)
  free_idx <- which(pt$free > 0L)

  sample_cov  <- lavaan::lavInspect(fit, "sampstat")$cov
  sample_nobs <- lavaan::lavInspect(fit, "nobs")

  imputations <- vector("list", m)

  for (i in seq_len(m)) {
    # level 1: fresh stochastic initial imputation
    x0_i <- validate_initializer_output(initializer(data), data, ov_names)

    # level 2: perturbed model-implied covariance from theta^(k)
    pt_i <- pt
    pt_i$start[free_idx] <- draws[i, ]
    pt_i$est[free_idx]   <- draws[i, ]

    # Suppress lavaan startup banner per iteration.  The model is not refit
    # (do.fit = FALSE): it is only instantiated to compute the model-implied
    # moments from the perturbed parameters, so the post-fit convergence
    # check must be disabled.
    fit_i <- suppressMessages(
      lavaan::lavaan(
        model      = pt_i,
        sample.cov = sample_cov,
        sample.nobs = sample_nobs,
        do.fit     = FALSE,
        se         = "none",
        test       = "none",
        check.post = FALSE
      )
    )

    sigma_i <- as.matrix(lavaan::fitted(fit_i)$cov)
    colnames(sigma_i) <- ov_names
    rownames(sigma_i) <- ov_names

    imputations[[i]] <- prism_project(
      data               = data,
      ov_names           = ov_names,
      sigma_target       = sigma_i,
      initial_imputation = x0_i,
      ...
    )
    attr(imputations[[i]], "imputation") <- i
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
