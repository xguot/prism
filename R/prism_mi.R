#' Proper Multiple Imputation via Parameter Perturbation
#'
#' @description Generate \code{m} completed datasets by drawing perturbed
#'   model-implied covariance matrices from the asymptotic sampling
#'   distribution of the FIML parameter estimates, then projecting the
#'   initial imputation onto each perturbed manifold. This produces proper
#'   multiple imputations (Rubin 1987) suitable for pooling with Rubin's
#'   rules.
#'
#' @details
#' The pipeline:
#' \enumerate{
#'   \item Extract free parameter estimates \eqn{\hat{\theta}} and their
#'         asymptotic covariance matrix \eqn{\text{ACOV}(\hat{\theta})}
#'         from the fitted lavaan object.
#'   \item Draw \code{m} perturbed parameter vectors from
#'         \eqn{\mathcal{N}(\hat{\theta}, \text{ACOV}(\hat{\theta}))}.
#'   \item For each draw, reconstruct the model-implied covariance matrix
#'         \eqn{\Sigma_k} by instantiating a lavaan model with the perturbed
#'         parameters and computing implied moments.
#'   \item Run the Lagrangian projection engine against each \eqn{\Sigma_k},
#'         producing \code{m} distinct completed datasets.
#' }
#'
#' @param data A data frame containing missing values.
#' @param fit A fitted lavaan object (must be estimated with
#'   \code{missing = "fiml"}).
#' @param m An integer specifying the number of imputations. Defaults to 20.
#' @param initial_imputation A data frame or matrix of the same dimensions as
#'   the longitudinal subset of \code{data}, containing initial imputed
#'   values. If \code{NULL}, column-mean imputation is used as a fallback.
#' @param return_mids A logical value. If \code{TRUE}, the result is cast to
#'   a \code{\link[mice]{mids}} object for compatibility with the
#'   \code{mice} pooling workflow. Requires the \code{mice} package.
#'   Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to the internal projection engine,
#'   e.g. \code{lambda}, \code{learning_rate}, \code{tol}, \code{max_iter}.
#'
#' @return A list of \code{m} completed data frames (class
#'   \code{"prism_mi_list"}), or a \code{mids} object if
#'   \code{return_mids = TRUE}.
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
#' mi_list <- prism_mi(df, fit, m = 20)
#' }
prism_mi <- function(data, fit, m = 20, initial_imputation = NULL,
                        return_mids = FALSE, ...) {
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.", call. = FALSE)
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for prism_mi(). ",
         "Install it with install.packages('MASS').", call. = FALSE)
  }

  theta_hat <- lavaan::coef(fit)
  acov <- lavaan::vcov(fit)

  if (anyNA(theta_hat) || anyNA(acov)) {
    stop("Parameter estimates or asymptotic covariance matrix contain NAs. ",
         "The model may not be identified.", call. = FALSE)
  }

  draws <- tryCatch(
    MASS::mvrnorm(n = m, mu = theta_hat, Sigma = acov),
    error = function(e) {
      # Fall back to perturbing each parameter independently
      matrix(rnorm(m * length(theta_hat), mean = theta_hat,
                   sd = sqrt(pmax(diag(acov), 1e-8))),
             nrow = m, byrow = TRUE)
    }
  )

  pt <- lavaan::parTable(fit)
  free_idx <- which(pt$free > 0L)

  sample_cov  <- lavaan::lavInspect(fit, "sampstat")$cov
  sample_nobs <- lavaan::lavInspect(fit, "nobs")
  time_cols   <- lavaan::lavNames(fit, "ov")

  imputations <- vector("list", m)

  for (i in seq_len(m)) {
    pt_i <- pt
    pt_i$start[free_idx] <- draws[i, ]
    pt_i$est[free_idx]   <- draws[i, ]

    # Suppress lavaan startup banner per iteration
    fit_i <- suppressMessages(
      lavaan::lavaan(
        model      = pt_i,
        sample.cov = sample_cov,
        sample.nobs = sample_nobs,
        do.fit     = FALSE,
        se         = "none",
        test       = "none"
      )
    )

    sigma_i <- lavaan::lavTech(fit_i, "implied")$cov
    colnames(sigma_i) <- time_cols
    rownames(sigma_i) <- time_cols

    imputations[[i]] <- prism_project(
      data               = data,
      time_cols          = time_cols,
      sigma_target       = sigma_i,
      initial_imputation = initial_imputation,
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
    # Remove the list class attribute so as.mids can assign its own
    imps_plain <- imputations
    class(imps_plain) <- "list"
    return(mice::as.mids(imps_plain))
  }

  imputations
}
