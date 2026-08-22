#' @title FIML Covariance Projection (PRISM v2)
#'
#' @description Fit a latent growth model via FIML, extract the model-implied
#'   covariance matrix, and project an initial imputation onto that structural
#'   target under the PRISM v2 regularized objective. Fidelity to the initial
#'   imputation (nonparametric, local information) is traded off against
#'   covariance matching (parametric structure) through \code{lambda_sigma},
#'   while the column means of the initial imputation are preserved exactly.
#'   Only originally-missing cells are modified; observed data are held fixed.
#'
#' @param data A data frame containing missing values.
#' @param model A lavaan model syntax string (e.g. a growth model).
#' @param initial_imputation A data frame or matrix of the same dimensions as
#'   the longitudinal subset of `data`, containing initial imputed values.
#'   If `NULL` (default), column-mean imputation is used as a fallback.
#' @param lambda_sigma A non-negative numeric value giving the structural
#'   weight of the covariance-matching term relative to the fidelity term.
#'   Both loss terms are normalized (fidelity per missing cell, covariance
#'   discrepancy per matrix entry), so \code{lambda_sigma} is dimensionless
#'   and comparable across sample sizes, variable counts, and missingness
#'   rates. Defaults to 1; a sensitivity grid of 0.25-4 is recommended. `0`
#'   returns the initial imputation; larger values enforce the FIML-implied
#'   structure more strongly.
#' @param lr A numeric value for the initial gradient descent step size.
#'   Because both loss terms are normalized, the natural step scale is O(1);
#'   the Armijo backtracking line search adapts it automatically. Defaults
#'   to 1.
#' @param tol_kkT A numeric value for the stationarity tolerance on the
#'   projected-gradient (KKT residual) norm. Defaults to 1e-4.
#' @param tol_cov A numeric value for the feasibility tolerance on the
#'   covariance gap to the target. Defaults to 1e-6.
#' @param max_iter An integer specifying the maximum number of iterations
#'   for the gradient descent projection. Defaults to 2000.
#' @param lambda Deprecated; use `lambda_sigma`.
#' @param learning_rate Deprecated; use `lr`.
#' @param tol Deprecated; use `tol_cov`.
#'
#' @details
#' The completed data solve the regularized constrained projection
#' \deqn{\min_X \frac{1}{2N_{\text{mis}}}\|M \odot (X - X^{(0)})\|_F^2 +
#' \frac{\lambda_\Sigma}{2p^2}\|\Sigma(X) - \Sigma_\text{FIML}\|_F^2}
#' subject to the observed cells being fixed and the column means of
#' \eqn{X^{(0)}} preserved. Both terms are normalized — fidelity per missing
#' cell and covariance discrepancy per matrix entry — so
#' \eqn{\lambda_\Sigma} is a dimensionless trade-off parameter. Mean
#' preservation is enforced structurally: the gradient of the missing cells
#' is centered at zero within each column before every step, so column sums
#' of the imputed values cannot drift during optimization. Columns with a
#' single missing cell are frozen by this constraint.
#'
#' Convergence is certified by first-order (KKT) stationarity of the
#' constrained problem, not by the raw gradient norm. The
#' \code{prism_diagnostics} attribute on the returned data frame reports the
#' KKT residual \code{r_kkT}, the covariance feasibility gap \code{feas_gap},
#' and the per-column Lagrange multiplier estimates \code{nu} of the mean
#' constraints. The \code{status} field classifies the termination as
#' \code{"converged_feasible"} (stationary and within the covariance
#' tolerance), \code{"converged"} (stationary with a nonzero covariance gap;
#' the gap reflects the fidelity/structure trade-off, not target
#' infeasibility), or a non-stationary stop (\code{"stalled_line_search"} or
#' \code{"max_iter_reached"}).
#'
#' @return A data frame with FIML-consistent, covariance-projected imputed
#'   values. Only the originally-missing cells are modified. The
#'   \code{prism_diagnostics} attribute contains the optimization diagnostics.
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
#' result <- prism_fiml(df, model)
#' attr(result, "prism_diagnostics")
#' }
prism_fiml <- function(data, model, initial_imputation = NULL,
                       lambda_sigma = NULL, lr = 1,
                       tol_kkT = 1e-4, tol_cov = 1e-6,
                       max_iter = 2000,
                       lambda = NULL, learning_rate = NULL, tol = NULL) {

  # backward-compatible mapping of deprecated arguments
  if (!is.null(lambda)) {
    warning("Argument 'lambda' is deprecated; use 'lambda_sigma'.",
            call. = FALSE)
    lambda_sigma <- lambda
  }
  if (!is.null(learning_rate)) {
    warning("Argument 'learning_rate' is deprecated; use 'lr'.",
            call. = FALSE)
    lr <- learning_rate
  }
  if (!is.null(tol)) {
    warning("Argument 'tol' is deprecated; use 'tol_cov'.", call. = FALSE)
    tol_cov <- tol
  }

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required. ",
         "Please install it with install.packages('lavaan').",
         call. = FALSE)
  }

  # Pre-flight validation of the model columns, before lavaan is invoked,
  # so that degenerate inputs fail with clear, deterministic messages
  model_cols <- validate_model_columns(data, model)

  # Fit growth model with FIML
  fit <- tryCatch(
    lavaan::growth(model, data = data, missing = "fiml"),
    error = function(e) {
      stop("lavaan FIML estimation failed: ", conditionMessage(e),
           call. = FALSE)
    }
  )

  # Extract the model-implied covariance and run the shared projection engine
  prism_from_fit(
    data               = data,
    fit                = fit,
    initial_imputation = initial_imputation,
    lambda_sigma       = lambda_sigma,
    lr                 = lr,
    tol_kkT            = tol_kkT,
    tol_cov            = tol_cov,
    max_iter           = max_iter
  )
}


# Pre-flight validation of the observed variables referenced by the model.
#
# Fail fast with deterministic, user-readable errors before lavaan sees the
# data: non-numeric columns, 100%-missing columns, and columns absent from
# the data are all caught here.  Latent variables are not data columns and
# are skipped.
#' @keywords internal
validate_model_columns <- function(data, model) {
  pt <- lavaan::lavParseModelString(model)
  model_vars <- unique(c(pt$lhs, pt$rhs))
  model_vars <- model_vars[nzchar(model_vars) & model_vars != "1"]
  ov <- model_vars[model_vars %in% names(data)]
  if (length(ov) == 0) {
    stop("No observed model variables found in data.", call. = FALSE)
  }

  non_num <- ov[!vapply(data[, ov, drop = FALSE], is.numeric, logical(1))]
  if (length(non_num) > 0) {
    stop("All model columns must be strictly numeric. Non-numeric column(s): ",
         paste(non_num, collapse = ", "), ".", call. = FALSE)
  }

  na_frac <- colSums(is.na(data[, ov, drop = FALSE])) / nrow(data)
  if (any(na_frac == 1)) {
    stop("Column(s) ", paste(names(which(na_frac == 1)), collapse = ", "),
         " are 100% missing; target covariance cannot be estimated.",
         call. = FALSE)
  }

  ov
}
