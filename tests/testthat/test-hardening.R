# Regression tests for the hardening fixes:
#   (a) the engine's n_mis-scaled initial step converges on large-n problems
#   (b) prism_draw_parameters eigenvalue flooring for degenerate ACOVs
#   (c) prism_mi's behavior on a degenerate (non-converged) fit

hd_sem_outlier <- function(n = 100, seed = 42) {
  set.seed(seed)
  # two-factor SEM with a structural regression, heavy-tailed indicators
  eta1 <- stats::rnorm(n, 0, 1)
  eta2 <- 0.5 * eta1 + stats::rnorm(n, 0, sqrt(0.75))
  lam <- c(0.7, 0.75, 0.8)
  gen <- function(eta) {
    out <- vapply(lam, function(l) l * eta + stats::rt(n, 1.5),
                  numeric(n))
    out[sample.int(n * length(lam), round(n * length(lam) * 0.10))] <-
      stats::rnorm(round(n * length(lam) * 0.10), 0, 12)   # outliers
    out
  }
  y <- cbind(gen(eta1), gen(eta2))
  colnames(y) <- paste0("y", 1:6)
  as.data.frame(y)
}

hd_inject_mnar <- function(df, miss) {
  # outcome-dependent missingness: large y4 values are missing with high
  # probability, which drives the FIML fit into non-convergence territory
  for (j in seq_len(ncol(df))) {
    x <- df[[j]]
    p <- miss * rank(-abs(x)) / length(x)
    na_idx <- runif(length(x)) < p
    df[na_idx, j] <- NA
  }
  df
}

hd_model <- function() {
  "eta1 =~ y1 + y2 + y3
   eta2 =~ y4 + y5 + y6
   eta2 ~ eta1"
}

test_that("engine converges on large-n problems with the default step", {
  set.seed(17)
  n <- 5000
  X0 <- matrix(stats::rnorm(n * 4), n, 4)
  # introduce correlation so the target is not trivial
  X0[, 2] <- X0[, 2] + 0.8 * X0[, 1]
  X0[, 3] <- X0[, 3] + 0.6 * X0[, 1] + 0.7 * X0[, 2]
  X0[, 4] <- X0[, 4] + 0.5 * X0[, 1] + 0.6 * X0[, 2] + 0.7 * X0[, 3]
  mask <- matrix(0, n, 4)
  mask[sample.int(n * 4, round(0.30 * n * 4))] <- 1

  r <- prism:::constrain_covariance_v2(
    X0, mask, diag(4),
    lambda_sigma = 10, lr = 1, max_iter = 500,
    tol_kkT = 1e-4, tol_cov = 1e-6
  )
  expect_true(r$converged)
  expect_true(r$status %in% c("converged", "converged_feasible"))
  expect_lt(r$iterations, 500)
  # the mean anchor (X0 sums) is preserved exactly
  expect_equal(colSums(r$X_refined * mask), colSums(X0 * mask),
               tolerance = 1e-8)
})

test_that("prism_draw_parameters floors degenerate eigenvalues", {
  set.seed(5)
  theta <- c(1, 0.5, 0.7, 0.8)
  # singular ACOV: third direction has zero variance
  v <- matrix(c(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0.5), 4, 4, byrow = TRUE)
  acov <- v %*% diag(c(0.2, 0.1, 0, 0.05)) %*% t(v)
  draws <- prism:::prism_draw_parameters(theta, acov, 2000)
  expect_true(all(is.finite(draws)))
  # the null direction is regularized, not frozen: draws have variance
  # around the floor, never exactly zero
  expect_gt(stats::sd(draws[, 3]), 1e-5)

  # well-conditioned ACOV: the empirical draw covariance recovers it
  acov_ok <- matrix(c(0.2, 0.05, 0.05, 0.1), 2, 2)
  draws_ok <- prism:::prism_draw_parameters(c(0, 0), acov_ok, 20000)
  expect_equal(stats::cov(draws_ok), acov_ok, tolerance = 0.02)
})

test_that("prism_mi degrades clearly on a non-converged fit", {
  skip_if_not_installed("lavaan")
  df <- hd_inject_mnar(hd_sem_outlier(100, seed = 20260901), 0.30)
  fit <- lavaan::sem(hd_model(), data = df, missing = "fiml")

  init <- df
  for (j in seq_len(ncol(init))) {
    na_idx <- is.na(init[[j]])
    init[[j]][na_idx] <- mean(init[[j]], na.rm = TRUE)
  }

  # deterministic level-1 + unavailable ACOV: no uncertainty source left
  expect_error(
    prism_mi(df, fit, m = 2, initial_imputation = init,
             lambda_sigma = 10, target_means = TRUE),
    "no source of between-imputation uncertainty"
  )

  # stochastic level-1: degrades to level-1-only with a warning
  hot_deck <- function(data) {
    out <- data
    for (j in seq_len(ncol(out))) {
      x <- out[[j]]
      na_idx <- is.na(x)
      if (any(na_idx)) {
        out[[j]][na_idx] <- sample(x[!na_idx], sum(na_idx), replace = TRUE)
      }
    }
    out
  }
  set.seed(3)
  mi <- suppressWarnings(prism_mi(
    df, fit, m = 2, initializer = hot_deck,
    lambda_sigma = 10, target_means = TRUE, max_iter = 2000
  ))
  expect_s3_class(mi, "prism_mi_list")
  expect_length(mi, 2)
  expect_true(all(vapply(mi, function(x) !anyNA(x), logical(1))))
})
