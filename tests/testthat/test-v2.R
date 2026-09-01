# helper: complete data from a linear latent growth model
sim_lgm_complete <- function(n = 200, seed = 2026) {
  set.seed(seed)
  int  <- stats::rnorm(n, 0, 1)
  slp  <- stats::rnorm(n, 0.3, 0.4)
  data.frame(
    T0 = int + stats::rnorm(n, 0, 0.6),
    T1 = int + 1 * slp + stats::rnorm(n, 0, 0.6),
    T2 = int + 2 * slp + stats::rnorm(n, 0, 0.6),
    T3 = int + 3 * slp + stats::rnorm(n, 0, 0.6)
  )
}

# helper: inject MCAR missingness column by column
inject_missing <- function(df, n_mis) {
  for (j in seq_along(n_mis)) {
    if (n_mis[j] > 0) {
      df[sample.int(nrow(df), n_mis[j]), j] <- NA
    }
  }
  df
}

lgm_model <- function() {
  "i =~ 1*T0 + 1*T1 + 1*T2 + 1*T3
   s =~ 0*T0 + 1*T1 + 2*T2 + 3*T3"
}

col_mean_impute <- function(df) {
  out <- df
  for (j in seq_len(ncol(out))) {
    na_idx <- is.na(out[[j]])
    if (any(na_idx)) {
      out[[j]][na_idx] <- mean(out[[j]], na.rm = TRUE)
    }
  }
  out
}

time_cols <- paste0("T", 0:3)

test_that("v2 projection preserves column means exactly", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    lambda_sigma = 1, max_iter = 3000)

  expect_equal(colMeans(res[, time_cols]),
               colMeans(init[, time_cols]),
               tolerance = 1e-8)
})

test_that("observed values are frozen exactly", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    lambda_sigma = 1, max_iter = 3000)

  df_mat  <- as.matrix(df[, time_cols])
  res_mat <- as.matrix(res[, time_cols])
  expect_equal(res_mat[!is.na(df_mat)], df_mat[!is.na(df_mat)],
               tolerance = 1e-12)
})

test_that("lambda_sigma trades fidelity against structure", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  fit <- lavaan::growth(lgm_model(), data = df, missing = "fiml")
  target <- lavaan::lavInspect(fit, "cov.ov")

  res_small <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                          lambda_sigma = 0.5, max_iter = 3000)
  res_big <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                        lambda_sigma = 10, max_iter = 3000)

  init_mat <- as.matrix(init[, time_cols])
  drift_small <- mean((as.matrix(res_small[, time_cols]) - init_mat)^2)
  drift_big   <- mean((as.matrix(res_big[, time_cols]) - init_mat)^2)
  expect_lt(drift_small, drift_big)

  gap_small <- sqrt(sum((stats::cov(res_small[, time_cols]) - target)^2))
  gap_big   <- sqrt(sum((stats::cov(res_big[, time_cols]) - target)^2))
  expect_lt(gap_big, gap_small)
})

test_that("lambda_sigma = 0 returns the initial imputation", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    lambda_sigma = 0, max_iter = 1000)

  expect_equal(as.matrix(res[, time_cols]), as.matrix(init[, time_cols]),
               tolerance = 1e-6)
})

test_that("KKT diagnostics are attached and classify termination", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    lambda_sigma = 1, max_iter = 3000)

  d <- attr(res, "prism_diagnostics")
  expect_true(is.list(d))
  for (field in c("status", "converged", "iterations", "r_kkT",
                  "feas_gap", "objective", "fidelity", "cov_term",
                  "nu", "grad_spread", "max_abs_nu", "n_mis",
                  "tol_cov")) {
    expect_true(field %in% names(d))
  }
  expect_true(d$status %in% c(
    "converged_feasible", "converged",
    "stalled_line_search", "max_iter_reached"
  ))
  expect_true(d$converged)
  expect_lte(d$r_kkT, 1e-3)
  expect_length(d$nu, 4)
})

test_that("columns with a single missing cell are frozen", {
  df <- inject_missing(sim_lgm_complete(), c(0, 1, 40, 50))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    lambda_sigma = 1, max_iter = 2000)

  missing_t1 <- which(is.na(df$T1))
  expect_equal(res$T1[missing_t1], init$T1[missing_t1], tolerance = 1e-10)
})

test_that("engine preserves missing-cell column sums and decreases objective", {
  set.seed(7)
  X0 <- matrix(stats::rnorm(240), 60, 4)
  mask <- matrix(0, 60, 4)
  mask[sample.int(240, 60)] <- 1
  target <- diag(4)

  r <- prism:::constrain_covariance_v2(
    X0, mask, target,
    lambda_sigma = 1, lr = 0.1, max_iter = 500,
    tol_kkT = 1e-6, tol_cov = 1e-6
  )

  expect_equal(colSums(r$X_refined * mask), colSums(X0 * mask),
               tolerance = 1e-10)

  # normalized objective: fcov = 0.5 * ||R||^2 / p^2, fidelity 0 at X0
  f_init <- (0.5 / 16) * sum((stats::cov(X0) - target)^2)
  expect_lte(r$objective, f_init + 1e-10)
  expect_true(r$status %in% c(
    "converged_feasible", "converged",
    "stalled_line_search", "max_iter_reached"
  ))
})

test_that("legacy v1 engine remains available", {
  set.seed(3)
  X0 <- matrix(stats::rnorm(90), 30, 3)
  mask <- matrix(0, 30, 3)
  mask[sample.int(90, 15)] <- 1
  r <- prism:::constrain_covariance(X0, mask, diag(3), 1.0, 0.001, 100, 1e-4)
  expect_true(is.matrix(r$X_refined))
  expect_true(r$converged || r$iterations >= 1)
})

test_that("prism_mi produces distinct two-level imputations reproducibly", {
  df <- inject_missing(sim_lgm_complete(80, seed = 99), c(0, 12, 16, 20))
  fit <- lavaan::growth(lgm_model(), data = df, missing = "fiml")

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

  set.seed(1)
  mi <- prism_mi(df, fit, m = 3, initializer = hot_deck, max_iter = 1500)

  expect_s3_class(mi, "prism_mi_list")
  expect_length(mi, 3)
  expect_true(all(vapply(mi, function(x) !anyNA(x), logical(1))))

  df_mat <- as.matrix(df[, time_cols])
  diff12 <- max(abs(as.matrix(mi[[1]][, time_cols]) -
                    as.matrix(mi[[2]][, time_cols])))
  diff13 <- max(abs(as.matrix(mi[[1]][, time_cols]) -
                    as.matrix(mi[[3]][, time_cols])))
  expect_gt(diff12, 0)
  expect_gt(diff13, 0)

  # missing-cell means vary across imputations (between-imputation variance)
  means <- vapply(mi, function(d) {
    x <- as.matrix(d[, time_cols])
    mean(x[is.na(df_mat)])
  }, numeric(1))
  expect_gt(stats::var(means), 0)

  # seed reproducibility
  set.seed(2)
  a <- prism_mi(df, fit, m = 2, initializer = hot_deck, max_iter = 1500)
  set.seed(2)
  b <- prism_mi(df, fit, m = 2, initializer = hot_deck, max_iter = 1500)
  expect_identical(a[[1]], b[[1]])
})

test_that("default bootstrap-forest initializer works", {
  skip_if_not_installed("missRanger")
  options(prism.num_trees = 50L)
  on.exit(options(prism.num_trees = NULL), add = TRUE)

  df <- inject_missing(sim_lgm_complete(60, seed = 5), c(0, 8, 12, 14))
  fit <- lavaan::growth(lgm_model(), data = df, missing = "fiml")

  set.seed(11)
  mi <- prism_mi(df, fit, m = 2, max_iter = 1000)
  expect_length(mi, 2)
  expect_gt(max(abs(as.matrix(mi[[1]][, time_cols]) -
                    as.matrix(mi[[2]][, time_cols]))), 0)
})

test_that("prism_mi converts to mids", {
  skip_if_not_installed("mice")
  df <- inject_missing(sim_lgm_complete(60, seed = 6), c(0, 8, 12, 14))
  fit <- lavaan::growth(lgm_model(), data = df, missing = "fiml")

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
  mids <- prism_mi(df, fit, m = 2, initializer = hot_deck,
                   return_mids = TRUE, max_iter = 800)
  expect_s3_class(mids, "mids")
  # complete(action = "long") returns the m imputations (original excluded)
  expect_equal(nrow(mice::complete(mids, action = "long")), nrow(df) * 2)
})

test_that("default lambda_sigma is 1 with normalized losses", {
  df <- inject_missing(sim_lgm_complete(100, seed = 8), c(0, 15, 20, 25))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    max_iter = 3000)
  d <- attr(res, "prism_diagnostics")
  expect_equal(d$lambda_sigma, 1)
  expect_equal(d$n_mis, sum(is.na(df[, time_cols])))
})

test_that("nonzero covariance gap at stationarity is not mislabeled", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                    lambda_sigma = 0.5, max_iter = 3000)
  d <- attr(res, "prism_diagnostics")
  expect_true(d$converged)
  # a finite lambda_sigma leaves a nonzero gap; this must be reported as a
  # regularized optimum, never as target infeasibility
  expect_true(d$status %in% c("converged", "converged_feasible"))
  expect_gt(d$feas_gap, 1e-4)
})

test_that("initial imputation observed cells are reset to raw data", {
  df <- inject_missing(sim_lgm_complete(), c(0, 30, 40, 50))
  init <- col_mean_impute(df)
  init[1, "T0"] <- init[1, "T0"] + 100   # corrupt an observed cell

  expect_warning(
    res <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                      lambda_sigma = 1, max_iter = 1000),
    "observed values"
  )
  expect_equal(res[1, "T0"], df[1, "T0"], tolerance = 1e-12)
  # anchors are computed after the reset, so the preserved means match the
  # corrected initial imputation
  corrected <- col_mean_impute(df)
  expect_equal(colMeans(res[, time_cols]),
               colMeans(corrected[, time_cols]), tolerance = 1e-8)
})

test_that("nearest_psd warns on a materially non-PSD input", {
  bad <- matrix(c(1, 1.5, 1.5, 1), 2, 2)  # eigenvalue -0.5, not numerical noise
  expect_warning(prism:::nearest_psd(bad), "PSD cone")
})

test_that("deprecated arguments warn and map to v2 parameters", {
  df <- inject_missing(sim_lgm_complete(100, seed = 8), c(0, 15, 20, 25))
  init <- col_mean_impute(df)

  expect_warning(
    res_old <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                          lambda = 0.5, learning_rate = 0.05, tol = 1e-3),
    "deprecated"
  )
  res_new <- prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
                        lambda_sigma = 0.5, lr = 0.05, tol_cov = 1e-3)
  expect_equal(res_old, res_new, tolerance = 1e-12)
})

test_that("invalid arguments are rejected", {
  df <- inject_missing(sim_lgm_complete(100, seed = 9), c(0, 15, 20, 25))
  init <- col_mean_impute(df)
  expect_error(
    prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init,
               lambda_sigma = -1),
    "non-negative"
  )
  expect_error(
    prism_fiml(df, lgm_model(), target_means = FALSE, initial_imputation = init, lr = 0),
    "positive"
  )
  fit <- lavaan::growth(lgm_model(), data = df, missing = "fiml")
  expect_error(prism_mi(df, fit, m = 1), "at least 2")
  expect_error(
    prism_mi(df, fit, m = 2, initializer = 1),
    "function"
  )
})
