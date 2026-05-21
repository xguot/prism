#' Smriti Multiple Imputation Wrapper
#'
#' This function generates multiple imputed datasets to capture the uncertainty
#' inherent in missing data imputation. It uses bootstrapping to create
#' distinct datasets, applies the Smriti Lagrangian refinement to each,
#' and returns a list of completed datasets suitable for pooling.
#'
#' @details
#' **Approximate (bootstrap) multiple imputation.** Proper multiple imputation
#' (Rubin 1987) draws imputations from the Bayesian posterior predictive
#' distribution of the missing data given the observed data and the imputation
#' model. This function instead bootstraps the rows and then deterministically
#' applies [smriti_impute()] to each bootstrap sample. The between-imputation
#' variance therefore captures sampling variability but *not* the full
#' imputation-model uncertainty. The resulting standard errors may be
#' moderately anti-conservative. Users needing valid Rubin's-rules pooling
#' should treat the output as approximate MI and consider adding a stochastic
#' residual-draw step, or use a fully Bayesian engine such as `mice`.
#'
#' Pooling the returned list can be performed manually via Rubin's rules:
#'
#' ```r
#' mi_list <- smriti_mi(df, time_cols = 2:5, m = 10)
#' # Within-imputation estimates
#' estimates <- lapply(mi_list, function(d) coef(lm(y ~ x, data = d)))
#' # Pool with Rubin's rules (point estimate = mean, SE = sqrt(W + (1+1/m)*B))
#' ```
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param m An integer specifying the number of imputations. Defaults to 5.
#' @param initial_imputation A matrix or data frame of the same dimensions as
#'   `data[, time_cols]`, containing initial imputed values. If `NULL`, the
#'   `smriti_impute` function's internal fallback (column-mean) is used.
#' @param ... Additional arguments passed to [smriti_impute()], e.g.
#'   `lambda`, `learning_rate`, `max_iter`, `robust`.
#'
#' @return A list of `m` completed data frames, each with an `"imputation"`
#'   attribute giving its index (1..m). The list has class `"smriti_mi_list"`.
#' @export
#'
#' @seealso [smriti_impute()] for the single-imputation engine.
#'
#' @examples
#' \dontrun{
#' # Assuming 'df' has missing longitudinal data in cols 2:5
#' mi_list <- smriti_mi(df, time_cols = 2:5, m = 5)
#' # Pool estimates manually via Rubin's rules
#' }
smriti_mi <- function(data, time_cols, m = 5, initial_imputation = NULL, ...) {
  n <- nrow(data)
  imputations <- vector("list", m)

  i <- 1
  attempts <- 0
  max_attempts <- m * 10

  while (i <= m) {
    attempts <- attempts + 1
    if (attempts > max_attempts) {
      stop(sprintf(
        "Failed to generate %d valid bootstrap samples after %d attempts. ",
        m, max_attempts
      ), "The data may be too sparse to support pairwise covariance ",
      "estimation across all time columns.")
    }

    # Bootstrap sample the rows to introduce variation for multiple imputation
    boot_idx <- sample(seq_len(n), n, replace = TRUE)
    boot_data <- data[boot_idx, , drop = FALSE]

    # Structural Safeguard: Verify that the bootstrap sample allows for
    # valid pairwise covariance/correlation estimation.
    cor_test <- suppressWarnings(
      stats::cor(boot_data[, time_cols], use = "pairwise.complete.obs")
    )
    if (anyNA(cor_test)) {
      next # Discard draw — a column pair has zero pairwise-complete cases
    }

    # Additional safeguard: check that the pairwise correlation matrix is
    # positive semidefinite (pairwise deletion can break this).
    cor_eig <- tryCatch(
      eigen(cor_test, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NULL
    )
    if (is.null(cor_eig) || min(cor_eig) < -1e-12) {
      next # Non-PSD correlation — discard and retry
    }

    # Handle initial imputation for the bootstrap sample
    boot_init <- NULL
    if (!is.null(initial_imputation)) {
      boot_init <- as.matrix(initial_imputation)[boot_idx, , drop = FALSE]
    }

    # Run the core smriti imputation on the bootstrapped data.
    # Wrap in tryCatch so that failures in smriti_impute (e.g. singular
    # target covariance) trigger a retry rather than aborting the MI loop.
    result <- tryCatch(
      smriti_impute(
        data              = boot_data,
        time_cols         = time_cols,
        initial_imputation = boot_init,
        ...
      ),
      error = function(e) NULL
    )
    if (is.null(result)) next

    # Tag the dataset with its imputation number
    attr(result, "imputation") <- i
    imputations[[i]] <- result
    i <- i + 1
  }

  class(imputations) <- c("smriti_mi_list", "list")
  return(imputations)
}
