# Mean targeting: the engine itself pins the missing-cell column sums so the
# returned matrix is the certified constrained optimum (joint mu + Sigma
# projection). Helpers are self-contained to keep this file independent of
# sibling test files.

mt_sim_lgm <- function(n = 200, seed = 2026) {
  set.seed(seed)
  int <- stats::rnorm(n, 0, 1)
  slp <- stats::rnorm(n, 0.3, 0.4)
  data.frame(
    T0 = int + stats::rnorm(n, 0, 0.6),
    T1 = int + 1 * slp + stats::rnorm(n, 0, 0.6),
    T2 = int + 2 * slp + stats::rnorm(n, 0, 0.6),
    T3 = int + 3 * slp + stats::rnorm(n, 0, 0.6)
  )
}

mt_inject <- function(df, n_mis) {
  for (j in seq_along(n_mis)) {
    if (n_mis[j] > 0) {
      df[sample.int(nrow(df), n_mis[j]), j] <- NA
    }
  }
  df
}

mt_lgm_model <- function() {
  "i =~ 1*T0 + 1*T1 + 1*T2 + 1*T3
   s =~ 0*T0 + 1*T1 + 2*T2 + 3*T3"
}

mt_col_mean <- function(df) {
  out <- df
  for (j in seq_len(ncol(out))) {
    na_idx <- is.na(out[[j]])
    if (any(na_idx)) {
      out[[j]][na_idx] <- mean(out[[j]], na.rm = TRUE)
    }
  }
  out
}

mt_time_cols <- paste0("T", 0:3)

test_that("target_means = TRUE matches model-implied means in columns with missing cells", {
  df <- mt_inject(mt_sim_lgm(200, seed = 4), c(0, 20, 25, 30))
  init <- mt_col_mean(df)
  fit <- lavaan::growth(mt_lgm_model(), data = df, missing = "fiml")
  mu_t <- lavaan::lavInspect(fit, "mean.ov")

  res <- prism_fiml(df, mt_lgm_model(), target_means = TRUE,
                    initial_imputation = init, lambda_sigma = 1,
                    max_iter = 3000)
  has_mis <- colSums(is.na(df)) > 0
  expect_equal(unname(colMeans(res)[has_mis]),
               unname(mu_t[has_mis]), tolerance = 1e-8)

  d <- attr(res, "prism_diagnostics")
  expect_equal(d$mean_gap, 0, tolerance = 1e-8)
  expect_equal(d$mu_target, unname(mu_t), tolerance = 1e-12)
})

test_that("a column with a single missing cell is pinned to the mean target", {
  df <- mt_inject(mt_sim_lgm(120, seed = 3), c(0, 1, 20, 25))
  init <- mt_col_mean(df)
  fit <- lavaan::growth(mt_lgm_model(), data = df, missing = "fiml")
  mu_t <- lavaan::lavInspect(fit, "mean.ov")

  res <- prism_fiml(df, mt_lgm_model(), target_means = TRUE,
                    initial_imputation = init, lambda_sigma = 1,
                    max_iter = 3000)

  # the single missing cell of T1 is fully determined by the mean anchor:
  # its value must be n * mu - observed sum
  mis_t1 <- which(is.na(df$T1))
  expect_equal(res$T1[mis_t1],
               unname(120 * mu_t["T1"] - sum(df$T1, na.rm = TRUE)),
               tolerance = 1e-6)
})

test_that("C++ engine holds s_j = n mu_j - obs_sum_j exactly", {
  set.seed(11)
  X0 <- matrix(stats::rnorm(300), 100, 3)
  mask <- matrix(0, 100, 3)
  mask[sample.int(300, 60)] <- 1
  target <- diag(3)
  mu_scaled <- c(0.3, -0.2, 0.5)

  r <- prism:::constrain_covariance_v2(
    X0, mask, target,
    lambda_sigma = 1, lr = 0.1, max_iter = 2000,
    tol_kkT = 1e-8, tol_cov = 1e-8,
    mu_scaled_target = mu_scaled
  )

  obs_sum <- colSums((1 - mask) * X0)
  anchor <- 100 * mu_scaled - obs_sum
  expect_equal(colSums(r$X_refined * mask), anchor, tolerance = 1e-8)
  expect_equal(colMeans(r$X_refined), mu_scaled, tolerance = 1e-8)
})

test_that("C++ engine falls back to the X0 anchor for NULL and NA entries", {
  set.seed(13)
  X0 <- matrix(stats::rnorm(240), 80, 3)
  mask <- matrix(0, 80, 3)
  mask[sample.int(240, 50)] <- 1
  target <- diag(3)

  # NULL target: legacy behaviour (column means of X0 preserved)
  r_null <- prism:::constrain_covariance_v2(
    X0, mask, target,
    lambda_sigma = 1, lr = 0.1, max_iter = 1000,
    tol_kkT = 1e-8, tol_cov = 1e-8
  )
  expect_equal(colSums(r_null$X_refined * mask), colSums(X0 * mask),
               tolerance = 1e-8)

  # per-entry NA keeps the X0 anchor for that column only
  r_mix <- prism:::constrain_covariance_v2(
    X0, mask, target,
    lambda_sigma = 1, lr = 0.1, max_iter = 1000,
    tol_kkT = 1e-8, tol_cov = 1e-8,
    mu_scaled_target = c(0.4, NA, NA)
  )
  expect_equal(sum(r_mix$X_refined[mask[, 1] == 1, 1]),
               80 * 0.4 - sum(X0[mask[, 1] == 0, 1]), tolerance = 1e-8)
  expect_equal(colSums(r_mix$X_refined[, 2:3] * mask[, 2:3]),
               colSums(X0[, 2:3] * mask[, 2:3]), tolerance = 1e-8)
})

test_that("target_means = FALSE still preserves X0 means exactly", {
  df <- mt_inject(mt_sim_lgm(200, seed = 8), c(0, 30, 40, 50))
  init <- mt_col_mean(df)
  res <- prism_fiml(df, mt_lgm_model(), target_means = FALSE,
                    initial_imputation = init, lambda_sigma = 1,
                    max_iter = 3000)
  expect_equal(colMeans(res[, mt_time_cols]),
               colMeans(init[, mt_time_cols]), tolerance = 1e-8)
  d <- attr(res, "prism_diagnostics")
  expect_true(is.na(d$mean_gap))
  expect_null(d$mu_target)
})

test_that("a custom numeric mean target is honored", {
  df <- mt_inject(mt_sim_lgm(150, seed = 7), c(0, 15, 20, 25))
  init <- mt_col_mean(df)
  custom <- c(0.5, 0.6, 0.7, 0.8)
  res <- prism_fiml(df, mt_lgm_model(), target_means = custom,
                    initial_imputation = init, max_iter = 3000)
  has_mis <- colSums(is.na(df)) > 0
  expect_equal(unname(colMeans(res)[has_mis]),
               unname(custom[has_mis]), tolerance = 1e-8)
})

test_that("target_means argument validation errors", {
  df <- mt_inject(mt_sim_lgm(100, seed = 6), c(0, 15, 20, 25))
  expect_error(
    prism_fiml(df, mt_lgm_model(), target_means = c(0, 1)),
    "number of observed variables"
  )
  expect_error(
    prism_fiml(df, mt_lgm_model(), target_means = c(0, 1, 2, NA)),
    "finite"
  )
  expect_error(
    prism_fiml(df, mt_lgm_model(), target_means = "yes"),
    "TRUE, FALSE, or a numeric vector"
  )
  expect_error(
    prism_fiml(df, mt_lgm_model(), target_means = c(TRUE, FALSE)),
    "TRUE, FALSE, or a numeric vector"
  )
})

test_that("diagnostics describe the returned matrix", {
  df <- mt_inject(mt_sim_lgm(200, seed = 5), c(0, 20, 25, 30))
  init <- mt_col_mean(df)
  fit <- lavaan::growth(mt_lgm_model(), data = df, missing = "fiml")
  sigma_t <- lavaan::lavInspect(fit, "cov.ov")
  mu_t <- lavaan::lavInspect(fit, "mean.ov")

  res <- prism:::prism_project(
    data = df, ov_names = mt_time_cols, sigma_target = sigma_t,
    mu_target = mu_t, initial_imputation = init,
    lambda_sigma = 1, max_iter = 3000
  )
  d <- attr(res, "prism_diagnostics")

  # rebuild the engine's scaled space from the initial imputation
  init_mat <- as.matrix(init)
  cmeans <- colMeans(init_mat)
  csds <- apply(init_mat, 2, stats::sd)
  X0_scaled <- scale(init_mat, center = cmeans, scale = csds)
  res_scaled <- scale(as.matrix(res), center = cmeans, scale = csds)

  # observed cells must round-trip untouched in the scaled space
  mask <- is.na(as.matrix(df))
  expect_equal(res_scaled[mask == 0], X0_scaled[mask == 0],
               tolerance = 1e-6)

  # ||S(returned) - Sigma_T|| recomputed in scaled space matches feas_gap
  Sigma_ret <- stats::cov(res_scaled)
  Sigma_scaled_target <- diag(1 / csds) %*% sigma_t %*% diag(1 / csds)
  gap_recomputed <- sqrt(sum((Sigma_ret - Sigma_scaled_target)^2))
  expect_equal(gap_recomputed, d$feas_gap, tolerance = 1e-6)

  # the mean target holds for the returned matrix itself
  has_mis <- colSums(mask) > 0
  expect_equal(unname(colMeans(as.matrix(res))[has_mis]),
               unname(mu_t[has_mis]), tolerance = 1e-8)
})
