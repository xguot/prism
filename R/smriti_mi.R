#' Smriti Multiple Imputation Wrapper
#'
#' This function generates multiple imputed datasets to capture the uncertainty
#' inherent in missing data imputation. It uses bootstrapping to create
#' distinct datasets, applies the Smriti Lagrangian refinement to each,
#' and returns a list of completed datasets suitable for pooling via Rubin's rules.
#'
#' @param data A data frame containing missing values.
#' @param time_cols A character vector or numeric vector specifying the
#'   longitudinal columns.
#' @param m An integer specifying the number of imputations. Defaults to 5.
#' @param initial_imputation A matrix or data frame of the same dimensions as 
#'   `data[, time_cols]`, containing initial imputed values. If `NULL`, the 
#'   `smriti_impute` function's internal fallback (column-mean) is used.
#' @param ... Additional arguments passed to [smriti_impute()].
#'
#' @return A list of `m` completed data frames.
#' @export
#' 
#' @examples
#' \dontrun{
#' # Assuming 'df' has missing longitudinal data in cols 2:5
#' mi_list <- smriti_mi(df, time_cols = 2:5, m = 5)
#' # Use mice::as.mids(mi_list) or similar to pool results
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
      stop(sprintf("Failed to generate %d valid bootstrap samples after %d attempts. ", m, max_attempts),
           "The data may be too sparse to support pairwise covariance estimation across all time columns.")
    }

    # Bootstrap sample the rows to introduce variation for multiple imputation
    boot_idx <- sample(seq_len(n), n, replace = TRUE)
    boot_data <- data[boot_idx, , drop = FALSE]
    
    # Structural Safeguard: Verify that the bootstrap sample allows for 
    # valid pairwise covariance/correlation estimation.
    cor_test <- suppressWarnings(stats::cor(boot_data[, time_cols], use = "pairwise.complete.obs"))
    if (anyNA(cor_test)) {
      next # Discard draw and try again
    }

    # Handle initial imputation for the bootstrap sample
    boot_init <- NULL
    if (!is.null(initial_imputation)) {
      boot_init <- as.matrix(initial_imputation)[boot_idx, , drop = FALSE]
    }

    # Run the core smriti imputation on the bootstrapped data
    imputations[[i]] <- smriti_impute(
      data = boot_data,
      time_cols = time_cols,
      initial_imputation = boot_init,
      ...
    )
    
    # Optional: tag the dataset with its imputation number
    attr(imputations[[i]], "imputation") <- i
    i <- i + 1
  }

  class(imputations) <- c("smriti_mi_list", "list")
  return(imputations)
}
