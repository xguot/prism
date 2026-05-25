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
