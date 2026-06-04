#' @title FIML-Smriti Refinement Wrapper
#'
#' @description This function serves as an agnostic refinement layer, combining
#'   Full Information Maximum Likelihood (FIML) with the
#'   structural covariance preservation of [smriti_impute()].
#'
#' @param data A data frame containing missing values.
#' @param model A lavaan model syntax string.
#' @param initial_imputation A data frame or matrix of the same dimensions as the
#'   subset of `data` defined by the model, containing initial imputed values.
#'   If `NULL` (default), `smriti_impute` handles initial imputation.
#' @param lambda A numeric value specifying the per-observation penalty weight
#'   for covariance matching. Defaults to 1.0.
#'
#' @return A data frame with FIML-consistent, covariance-projected imputed values.
#' @export
smriti_fiml <- function(data, model, initial_imputation = NULL, lambda = 1.0) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for this wrapper. ",
         "Please install it with install.packages('lavaan').",
         call. = FALSE)
  }

  # ── Fit growth model with FIML ──────────────────────────────────────────
  fit <- tryCatch(
    lavaan::growth(model, data = data, missing = "fiml"),
    error = function(e) {
      stop("lavaan FIML estimation failed: ", e$message, call. = FALSE)
    }
  )

  # ── Extract model-implied covariance as the structural target ────────────
  sigma_target <- tryCatch(
    lavaan::lavInspect(fit, "cov.ov"),
    error = function(e) {
      stop("Failed to extract model-implied covariance from lavaan fit: ",
           e$message, call. = FALSE)
    }
  )

  time_cols <- colnames(sigma_target)

  # ── Delegate to smriti_impute with FIML-derived target ───────────────────
  smriti_impute(
    data               = data,
    time_cols          = time_cols,
    initial_imputation = initial_imputation,
    lambda             = lambda,
    custom_target      = sigma_target,
    robust             = FALSE
  )
}
