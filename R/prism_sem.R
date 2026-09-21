#' @title SEM Covariance Projection (PRISM v2)
#'
#' @description Fit a structural equation model via FIML (or reuse an already
#'   fitted model), extract the model-implied covariance matrix, and project
#'   an initial imputation onto that structural target under the PRISM v2
#'   regularized objective. The engine, objective, mean anchoring, and KKT
#'   diagnostics are identical to \code{\link{prism_fiml}}; this entry point
#'   supports any single-group SEM estimated by \code{lavaan::sem()} —
#'   CFAs, structural regressions among latent variables, mediation models,
#'   and models with covariates. By default the completed column means target
#'   the FIML model-implied means. Only originally-missing cells are
#'   modified; observed data are held fixed.
#'
#' @param data A data frame containing missing values.
#' @param model Either lavaan model syntax (a single string) or a fitted
#'   lavaan object. When syntax is supplied it is fit with
#'   \code{lavaan::sem(model, data = data, missing = "fiml")}. When a fitted
#'   object is supplied it is used as-is and should have been estimated with
#'   \code{missing = "fiml"}.
#' @param initial_imputation A data frame or matrix of the same dimensions as
#'   the observed variables of the model, containing initial imputed values.
#'   If `NULL` (default), column-mean imputation is used as a fallback.
#' @param lambda_sigma A non-negative numeric value giving the structural
#'   weight of the covariance-matching term relative to the fidelity term.
#'   Both loss terms are normalized (fidelity per missing cell, covariance
#'   discrepancy per matrix entry), so \code{lambda_sigma} is dimensionless
#'   and comparable across sample sizes, variable counts, and missingness
#'   rates. Defaults to 1; a sensitivity grid of 0.25-4 is recommended. `0`
#'   returns the initial imputation; larger values enforce the
#'   model-implied structure more strongly.
#' @param lr A dimensionless multiplier of the natural initial step size.
#'   The engine starts each Armijo line search at \code{lr * N_mis}, the
#'   natural scale of the normalized objective (for the fidelity term alone
#'   this is the exact line minimum), and backtracks as needed. Defaults
#'   to 1.
#' @param tol_kkT A numeric value for the stationarity tolerance on the
#'   projected-gradient (KKT residual) norm. Defaults to 1e-4.
#' @param tol_cov A numeric value for the feasibility tolerance on the
#'   covariance gap to the target. Defaults to 1e-6.
#' @param max_iter An integer specifying the maximum number of iterations
#'   for the gradient descent projection. Defaults to 2000.
#' @param target_means A logical or numeric vector controlling the mean
#'   target of the projection. `TRUE` (default) targets the FIML
#'   model-implied means of the observed variables, so the engine solves a
#'   joint mean + covariance projection consistent under MAR. `FALSE`
#'   preserves the column means of the initial imputation exactly (legacy
#'   behaviour). A numeric vector of length p supplies a custom mean target.
#'   Columns without missing values keep their observed means in all modes.
#'
#' @details
#' The completed data solve the same regularized projection as
#' \code{\link{prism_fiml}}: fidelity to the initial imputation is traded off
#' against covariance matching through \code{lambda_sigma}, and each column
#' mean is fixed to its target — the FIML model-implied means by default
#' (\code{target_means = TRUE}) or the initial imputation's column means
#' when \code{target_means = FALSE}. Convergence is certified by the KKT
#' diagnostics in the \code{prism_diagnostics} attribute (see
#' \code{\link{prism_fiml}} for the field descriptions).
#'
#' Restrictions: the model must be single-group and single-level with
#' continuous observed variables. Ordinal, multi-group, and multilevel models
#' are rejected with an error.
#'
#' @return A data frame with the missing cells completed. Only the
#'   originally-missing cells are modified. The \code{prism_diagnostics}
#'   attribute contains the optimization diagnostics.
#' @export
#'
#' @examples
#' \dontrun{
#' library(lavaan)
#' model <- "
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#' "
#' df <- lavaan::PoliticalDemocracy
#' df$x1[1:10] <- NA
#' result <- prism_sem(df, model)
#' attr(result, "prism_diagnostics")
#' }
prism_sem <- function(data, model, initial_imputation = NULL,
                      lambda_sigma = NULL, lr = 1,
                      tol_kkT = 1e-4, tol_cov = 1e-6,
                      max_iter = 2000,
                      target_means = TRUE) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required. ",
         "Please install it with install.packages('lavaan').",
         call. = FALSE)
  }

  if (inherits(model, "lavaan")) {
    fit <- model
    validate_fitted_model(fit, data)
  } else if (is.character(model) && length(model) == 1L && nzchar(model)) {
    validate_model_columns(data, model)
    fit <- tryCatch(
      lavaan::sem(model, data = data, missing = "fiml"),
      error = function(e) {
        stop("lavaan FIML estimation failed: ", conditionMessage(e),
             call. = FALSE)
      }
    )
  } else {
    stop("'model' must be lavaan model syntax or a fitted lavaan object.",
         call. = FALSE)
  }

  prism_from_fit(
    data               = data,
    fit                = fit,
    initial_imputation = initial_imputation,
    lambda_sigma       = lambda_sigma,
    lr                 = lr,
    tol_kkT            = tol_kkT,
    tol_cov            = tol_cov,
    max_iter           = max_iter,
    target_means       = target_means
  )
}


# Validate a user-supplied fitted lavaan object against the supported scope:
# single-group, single-level, continuous observed variables present in data,
# and (warn-only) FIML estimation.
#' @keywords internal
validate_fitted_model <- function(fit, data) {
  if (lavaan::lavInspect(fit, "ngroups") != 1L) {
    stop("Multi-group SEM is not supported: project each group separately.",
         call. = FALSE)
  }
  if (lavaan::lavInspect(fit, "nlevels") != 1L) {
    stop("Multilevel SEM is not supported.", call. = FALSE)
  }
  ordered_vars <- lavaan::lavInspect(fit, "ordered")
  if (length(ordered_vars) > 0L) {
    stop("Ordinal observed variables are not supported: the continuous ",
         "covariance projection requires numeric outcomes. Affected: ",
         paste(ordered_vars, collapse = ", "), ".", call. = FALSE)
  }

  ov_names <- lavaan::lavNames(fit, "ov")
  missing_cols <- setdiff(ov_names, names(data))
  if (length(missing_cols) > 0L) {
    stop("Observed variables not found in data: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  non_num <- ov_names[!vapply(data[, ov_names, drop = FALSE],
                              is.numeric, logical(1))]
  if (length(non_num) > 0L) {
    stop("All observed variables must be strictly numeric. Non-numeric: ",
         paste(non_num, collapse = ", "), ".", call. = FALSE)
  }

  if (!isTRUE(identical(fit@Options$missing, "fiml"))) {
    warning("The supplied fit was not estimated with missing = \"fiml\"; ",
            "the structural target may be inconsistent with the ",
            "missing-data mechanism.", call. = FALSE)
  }
  invisible(NULL)
}


# Shared pipeline used by prism_fiml and prism_sem: extract the
# model-implied covariance (and optionally the means) of the observed
# variables from a fitted lavaan object and run the projection engine.
#' @keywords internal
prism_from_fit <- function(data, fit, initial_imputation,
                           lambda_sigma, lr, tol_kkT, tol_cov, max_iter,
                           target_means = TRUE) {
  sigma_target <- tryCatch(
    lavaan::lavInspect(fit, "cov.ov"),
    error = function(e) {
      stop("Failed to extract model-implied covariance from lavaan fit: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  ov_names <- colnames(sigma_target)

  # Resolve the mean target -----------------------------------------------
  # TRUE  -> FIML model-implied means (the default; consistent under MAR,
  #          removes the initialiser's first-moment bias)
  # FALSE -> mean preservation of the initial imputation (legacy behaviour)
  # numeric vector -> user-supplied custom mean target
  implied_means <- tryCatch(
    lavaan::lavInspect(fit, "mean.ov"),
    error = function(e) NULL
  )
  mu_target <- resolve_mean_target(target_means, implied_means, ov_names)

  prism_project(
    data               = data,
    ov_names           = ov_names,
    sigma_target       = sigma_target,
    mu_target          = mu_target,
    initial_imputation = initial_imputation,
    lambda_sigma       = lambda_sigma,
    lr                 = lr,
    tol_kkT            = tol_kkT,
    tol_cov            = tol_cov,
    max_iter           = max_iter
  )
}
