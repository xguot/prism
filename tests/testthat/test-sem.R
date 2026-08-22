# helper: two-factor SEM with a structural regression among latents
sim_sem_complete <- function(n = 300, seed = 2026) {
  set.seed(seed)
  eta1 <- stats::rnorm(n, 0, 1)
  eta2 <- 0.5 * eta1 + stats::rnorm(n, 0, sqrt(0.75))
  lam <- c(0.7, 0.75, 0.8)
  gen <- function(eta) {
    vapply(lam, function(l) l * eta + stats::rnorm(n, 0, sqrt(1 - l^2)),
           numeric(n))
  }
  y <- cbind(gen(eta1), gen(eta2))
  colnames(y) <- paste0("y", 1:6)
  as.data.frame(y)
}

sem_model <- function() {
  "eta1 =~ y1 + y2 + y3
   eta2 =~ y4 + y5 + y6
   eta2 ~ eta1"
}

inject_missing_sem <- function(df, n_mis) {
  for (j in seq_along(n_mis)) {
    if (n_mis[j] > 0) {
      df[sample.int(nrow(df), n_mis[j]), j] <- NA
    }
  }
  df
}

col_mean_impute_sem <- function(df) {
  out <- df
  for (j in seq_len(ncol(out))) {
    na_idx <- is.na(out[[j]])
    if (any(na_idx)) {
      out[[j]][na_idx] <- mean(out[[j]], na.rm = TRUE)
    }
  }
  out
}

hot_deck_sem <- function(data) {
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

test_that("prism_sem completes a CFA with structural paths", {
  df <- inject_missing_sem(sim_sem_complete(), c(40, 45, 50, 42, 38, 44))
  init <- col_mean_impute_sem(df)
  res <- prism_sem(df, sem_model(), initial_imputation = init,
                   lambda_sigma = 50, max_iter = 3000)

  expect_false(anyNA(res))
  expect_equal(dim(res), dim(df))

  df_mat  <- as.matrix(df)
  res_mat <- as.matrix(res)
  expect_equal(res_mat[!is.na(df_mat)], df_mat[!is.na(df_mat)],
               tolerance = 1e-12)
  expect_equal(colMeans(res), colMeans(init), tolerance = 1e-8)

  d <- attr(res, "prism_diagnostics")
  expect_true(d$status %in% c(
    "converged_feasible", "constrained_geometric_limit",
    "geometric_infeasible", "stalled_line_search", "max_iter_reached"
  ))
  expect_true(d$converged)
  expect_length(d$nu, 6)
})

test_that("prism_sem accepts a pre-fitted lavaan object", {
  df <- inject_missing_sem(sim_sem_complete(), c(40, 45, 50, 42, 38, 44))
  init <- col_mean_impute_sem(df)
  fit <- lavaan::sem(sem_model(), data = df, missing = "fiml")

  res_fit <- prism_sem(df, fit, initial_imputation = init,
                       lambda_sigma = 50, max_iter = 3000)
  res_syn <- prism_sem(df, sem_model(), initial_imputation = init,
                       lambda_sigma = 50, max_iter = 3000)
  expect_equal(res_fit, res_syn, tolerance = 1e-12)
})

test_that("prism_sem rejects ordinal models", {
  df <- inject_missing_sem(sim_sem_complete(200), c(25, 25, 25, 25, 25, 25))
  # discretize one indicator into four categories, then fit with WLSMV
  df$y1o <- as.ordered(cut(df$y1, breaks = stats::quantile(df$y1, c(0, .25, .5, .75, 1), na.rm = TRUE),
                           include.lowest = TRUE, labels = FALSE))
  model_ord <- "eta1 =~ y1o + y2 + y3
                eta2 =~ y4 + y5 + y6
                eta2 ~ eta1"
  fit_ord <- lavaan::sem(model_ord, data = df, ordered = "y1o")
  expect_error(prism_sem(df, fit_ord), "rdinal")
})

test_that("prism_sem rejects multi-group models", {
  df <- inject_missing_sem(sim_sem_complete(), c(30, 30, 30, 30, 30, 30))
  df$g <- rep(c("a", "b"), length.out = nrow(df))
  fit_mg <- lavaan::sem(sem_model(), data = df, group = "g",
                        missing = "fiml")
  expect_error(prism_sem(df, fit_mg), "ulti-group")
})

test_that("prism_sem warns when the supplied fit was not FIML-estimated", {
  df <- inject_missing_sem(sim_sem_complete(), c(30, 30, 30, 30, 30, 30))
  init <- col_mean_impute_sem(df)
  fit_lw <- lavaan::sem(sem_model(), data = df)  # listwise by default
  expect_warning(
    prism_sem(df, fit_lw, initial_imputation = init, max_iter = 1000),
    "fiml"
  )
})

test_that("prism_sem rejects invalid model arguments", {
  df <- inject_missing_sem(sim_sem_complete(), c(30, 30, 30, 30, 30, 30))
  expect_error(prism_sem(df, 42), "syntax or a fitted")
})

test_that("prism_mi works on SEM fits", {
  df <- inject_missing_sem(sim_sem_complete(120, seed = 7),
                           c(15, 18, 20, 16, 19, 21))
  fit <- lavaan::sem(sem_model(), data = df, missing = "fiml")

  set.seed(1)
  mi <- prism_mi(df, fit, m = 2, initializer = hot_deck_sem,
                 max_iter = 1000)
  expect_s3_class(mi, "prism_mi_list")
  expect_length(mi, 2)
  expect_true(all(vapply(mi, function(x) !anyNA(x), logical(1))))
  expect_gt(max(abs(as.matrix(mi[[1]]) - as.matrix(mi[[2]]))), 0)
})

test_that("stochastic_fiml_impute works on SEM fits", {
  df <- inject_missing_sem(sim_sem_complete(120, seed = 8),
                           c(15, 18, 20, 16, 19, 21))
  fit <- lavaan::sem(sem_model(), data = df, missing = "fiml")
  res <- stochastic_fiml_impute(df, fit)

  expect_false(anyNA(res))
  df_mat  <- as.matrix(df)
  res_mat <- as.matrix(res)
  expect_equal(res_mat[!is.na(df_mat)], df_mat[!is.na(df_mat)],
               tolerance = 1e-12)
})
