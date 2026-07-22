#' @title FIML Covariance Projection
#'
#' @description Fit a latent growth model via FIML, extract the model-implied
#'   covariance matrix, and project an initial imputation onto that structural
#'   manifold using Lagrangian-constrained gradient descent. Only
#'   originally-missing cells are modified; observed data is held fixed.
#'
#' @param data A data frame containing missing values.
#' @param model A lavaan model syntax string (e.g. a growth model).
#' @param initial_imputation A data frame or matrix of the same dimensions as
#'   the longitudinal subset of `data`, containing initial imputed values.
#'   If `NULL` (default), column-mean imputation is used as a fallback.
#' @param lambda A numeric value specifying the per-observation penalty weight
#'   for covariance matching. Defaults to 1.0.
#' @param learning_rate A numeric value for the gradient descent step size.
#'   Defaults to 0.001.
#' @param tol A numeric value for the convergence tolerance (Frobenius norm).
#'   Defaults to 1e-6.
#' @param max_iter An integer specifying the maximum number of iterations
#'   for the gradient descent projection. Defaults to 2000.
#'
#' @return A data frame with FIML-consistent, covariance-projected imputed
#'   values. Only the originally-missing cells are modified.
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
#' prism_fiml(df, model)
#' }
prism_fiml <- function(data, model, initial_imputation = NULL,
                         lambda = 1.0, learning_rate = 0.001,
                         tol = 1e-6, max_iter = 2000) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required. ",
         "Please install it with install.packages('lavaan').",
         call. = FALSE)
  }

  # Fit growth model with FIML
  fit <- tryCatch(
    lavaan::growth(model, data = data, missing = "fiml"),
    error = function(e) {
      stop("lavaan FIML estimation failed: ", e$message, call. = FALSE)
    }
  )

  # Extract model-implied covariance as the structural target
  sigma_target <- tryCatch(
    lavaan::lavInspect(fit, "cov.ov"),
    error = function(e) {
      stop("Failed to extract model-implied covariance from lavaan fit: ",
           e$message, call. = FALSE)
    }
  )

  time_cols <- colnames(sigma_target)

  # Delegate to the internal projection engine
  prism_project(
    data               = data,
    time_cols          = time_cols,
    sigma_target       = sigma_target,
    initial_imputation = initial_imputation,
    lambda             = lambda,
    learning_rate      = learning_rate,
    tol                = tol,
    max_iter           = max_iter
  )
}
