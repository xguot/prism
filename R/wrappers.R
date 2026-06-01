#' @title missForest-Smriti Refinement Wrapper
#'
#' @description This function serves as an agnostic refinement layer, combining
#'   the predictive power of Random Forests via `missForest` with the
#'   structural covariance preservation of [smriti_impute()].
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param robust A logical value. If TRUE (default), uses robust covariance estimation.
#' @param ... Additional arguments passed directly to [missForest::missForest()].
#'
#' @return A data frame with covariance-projected imputed values.
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4))
#' smriti_forest(df, time_cols = 1:2)
#' }
smriti_forest <- function(data, time_cols, robust = TRUE, ...) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("Package 'missForest' is required for this wrapper. ",
         "Please install it with install.packages('missForest').",
         call. = FALSE)
  }

  # Extract subset for imputation
  x_subset <- data[, time_cols, drop = FALSE]

  # Initial imputation via missForest
  forest_res <- missForest::missForest(x_subset, ...)
  ximp <- forest_res$ximp

  # Refine via smriti
  smriti_impute(data = data, time_cols = time_cols, initial_imputation = ximp, robust = robust)
}

#' @title missRanger-Smriti Refinement Wrapper
#'
#' @description This function serves as an agnostic refinement layer, combining
#'   the fast predictive imputation of `missRanger` with the structural
#'   covariance preservation of [smriti_impute()].
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param robust A logical value. If TRUE (default), uses robust covariance estimation.
#' @param ... Additional arguments passed directly to [missRanger::missRanger()].
#'
#' @return A data frame with covariance-projected imputed values.
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4))
#' smriti_ranger(df, time_cols = 1:2)
#' }
smriti_ranger <- function(data, time_cols, robust = TRUE, ...) {
  if (!requireNamespace("missRanger", quietly = TRUE)) {
    stop("Package 'missRanger' is required for this wrapper. ",
         "Please install it with install.packages('missRanger').",
         call. = FALSE)
  }

  # Extract subset for imputation
  x_subset <- data[, time_cols, drop = FALSE]

  # Initial imputation via missRanger
  ximp <- missRanger::missRanger(x_subset, ...)

  # Refine via smriti
  smriti_impute(data = data, time_cols = time_cols, initial_imputation = ximp, robust = robust)
}

#' @title FIML-Smriti Refinement Wrapper
#'
#' @description Bridges the gap between Full Information Maximum Likelihood
#'   (FIML) estimation and complete-data analysis.  A latent growth curve
#'   model is fitted to the incomplete data via [lavaan::growth()] with
#'   `missing = "fiml"`, the model-implied covariance matrix is extracted,
#'   and [smriti_impute()] uses it as a `custom_target` to produce a
#'   completed dataset whose covariance structure matches the FIML estimate.
#'   This enables FIML-quality covariance in any downstream model that
#'   requires complete data (random forests, PCA, clustering, etc.).
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param model A lavaan model syntax string.  Defaults to a linear growth
#'   model with fixed time scores `0, 1, ..., p-1`.  If `NULL`, a linear
#'   growth model is auto-generated from the number of `time_cols`.
#' @param ... Additional arguments passed to [smriti_impute()], e.g.
#'   `lambda`, `learning_rate`, `max_iter`.
#'
#' @return A data frame with FIML-consistent, covariance-projected imputed values.
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4), T3 = c(1, 2, NA, 4))
#' smriti_fiml(df, time_cols = c("T1", "T2", "T3"))
#' }
smriti_fiml <- function(data, time_cols, model = NULL, ...) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for this wrapper. ",
         "Please install it with install.packages('lavaan').",
         call. = FALSE)
  }

  time_names <- colnames(data[, time_cols, drop = FALSE])
  p <- length(time_names)

  if (p < 2) {
    stop("At least two time columns are required for a growth model.",
         call. = FALSE)
  }

  # ── Auto-generate linear growth model if none supplied ───────────────────
  if (is.null(model)) {
    rhs <- paste(time_names, collapse = " + ")
    time_scores <- seq_len(p) - 1
    i_loadings <- paste0("1*", time_names, collapse = " + ")
    s_loadings <- paste0(time_scores, "*", time_names, collapse = " + ")
    model <- paste0(
      "i =~ ", i_loadings, "\n",
      "s =~ ", s_loadings, "\n",
      "i ~~ s\n",
      "i ~~ i\n",
      "s ~~ s\n"
    )
  }

  # ── Fit growth model with FIML ──────────────────────────────────────────
  fit <- tryCatch(
    lavaan::growth(model, data = data, missing = "fiml",
                   fixed.x = FALSE, auto.var = TRUE),
    error = function(e) {
      stop("lavaan FIML estimation failed: ", e$message, call. = FALSE)
    }
  )

  if (!lavaan::lavInspect(fit, "converged")) {
    warning("lavaan model did not converge.  ",
            "The model-implied covariance may be unreliable.")
  }

  # ── Extract model-implied covariance as the structural target ────────────
  sigma_target <- tryCatch(
    lavaan::lavInspect(fit, "cov.ov"),
    error = function(e) {
      stop("Failed to extract model-implied covariance from lavaan fit: ",
           e$message, call. = FALSE)
    }
  )

  # ── Delegate to smriti_impute with FIML-derived target ───────────────────
  smriti_impute(
    data          = data,
    time_cols     = time_cols,
    custom_target = sigma_target,
    robust        = FALSE,
    ...
  )
}

#' @title mice-Smriti Refinement Wrapper
#'
#' @description This function serves as an agnostic refinement layer, combining
#'   multivariate imputation by chained equations (MICE) with the
#'   structural covariance preservation of [smriti_impute()].
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param robust A logical value. If TRUE (default), uses robust covariance estimation.
#' @param ... Additional arguments passed directly to [mice::mice()].
#'
#' @return A data frame with covariance-projected imputed values.
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(T1 = c(1, NA, 3, 4), T2 = c(NA, 2, 3, 4))
#' smriti_mice(df, time_cols = 1:2)
#' }
smriti_mice <- function(data, time_cols, robust = TRUE, ...) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for this wrapper. ",
         "Please install it with install.packages('mice').",
         call. = FALSE)
  }

  # Extract subset for imputation
  x_subset <- data[, time_cols, drop = FALSE]

  # Initial imputation via mice
  imp <- mice::mice(x_subset, ...)

  # Extract first completed dataset
  ximp <- mice::complete(imp, 1)

  # Refine via smriti
  smriti_impute(data = data, time_cols = time_cols, initial_imputation = ximp, robust = robust)
}
