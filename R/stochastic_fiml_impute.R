#' Stochastic Regression Imputation via FIML Parameters
#'
#' @description Traditional single imputation baseline: for each incomplete
#'   row, compute the conditional expectation of missing variables given
#'   observed variables under the model-implied multivariate normal
#'   distribution, then add a stochastic residual draw to preserve variance.
#'   This is the classic method (Schafer 1997) — it preserves first and
#'   second moments asymptotically but does not target the exact sample
#'   covariance structure.
#'
#' @param data A data frame containing missing values.
#' @param fit A fitted lavaan object (estimated with
#'   \code{missing = "fiml"}).
#'
#' @return A single completed data frame with observed values untouched and
#'   missing values replaced by stochastic draws from the conditional
#'   predictive distribution.
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
#' completed <- stochastic_fiml_impute(df, fit)
#' }
stochastic_fiml_impute <- function(data, fit) {
  if (!inherits(fit, "lavaan")) {
    stop("'fit' must be a fitted lavaan object.", call. = FALSE)
  }

  vars <- lavaan::lavNames(fit, "ov")
  implied <- lavaan::lavTech(fit, "implied")
  mu_all <- as.vector(implied$mean)
  Sigma_all <- as.matrix(implied$cov)
  p <- length(vars)

  if (is.null(mu_all)) {
    mu_all <- rep(0, p)
  }

  x <- as.matrix(data[, vars])
  na_rows <- which(rowSums(is.na(x)) > 0L)

  for (i in na_rows) {
    x_i <- x[i, ]
    obs_idx <- which(!is.na(x_i))
    mis_idx <- which(is.na(x_i))

    if (length(obs_idx) == 0L) {
      x[i, mis_idx] <- tryCatch(
        MASS::mvrnorm(n = 1L, mu = mu_all[mis_idx],
                      Sigma = as.matrix(Sigma_all[mis_idx, mis_idx, drop = FALSE])),
        error = function(e) mu_all[mis_idx]
      )
      next
    }

    mu_o  <- mu_all[obs_idx]
    mu_m  <- mu_all[mis_idx]
    S_oo  <- as.matrix(Sigma_all[obs_idx, obs_idx, drop = FALSE])
    S_mm  <- as.matrix(Sigma_all[mis_idx, mis_idx, drop = FALSE])
    S_mo  <- as.matrix(Sigma_all[mis_idx, obs_idx, drop = FALSE])
    S_oo_inv <- tryCatch(solve(S_oo), error = function(e) {
      tryCatch(MASS::ginv(S_oo), error = function(e2) {
        solve(S_oo + diag(1e-8, nrow(S_oo)))
      })
    })

    mu_cond <- mu_m + S_mo %*% S_oo_inv %*% (x_i[obs_idx] - mu_o)
    S_cond  <- S_mm - S_mo %*% S_oo_inv %*% t(S_mo)
    S_cond  <- (S_cond + t(S_cond)) / 2

    x[i, mis_idx] <- tryCatch(
      MASS::mvrnorm(n = 1L, mu = as.vector(mu_cond), Sigma = S_cond),
      error = function(e) as.vector(mu_cond)
    )
  }

  data[, vars] <- x
  data
}
