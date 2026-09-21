# Smoke test for the GCM + SEM simulation workers.
#
# Run from the simulations/ directory with the package library on the path:
#   R_LIBS=../libs Rscript smoke_test.R
#
# Sources the simulation definitions (everything above the SLURM dispatch
# block), then runs one replication of each model worker and checks the
# output schema and the key PRISM invariants.

.libPaths(c(normalizePath(file.path("..", "libs")), .libPaths()))

lines <- readLines("simulation.R")
cut <- grep("SLURM Array Dispatch", lines)[1] - 1L
eval(parse(text = paste(lines[seq_len(cut)], collapse = "\n")), envir = .GlobalEnv)

stopifnot(two_level_mi == FALSE)   # deterministic level-1 path by default

check_schema <- function(res, model, label) {
  cat(sprintf("%s worker returned %d rows x %d cols\n",
              label, nrow(res), ncol(res)))
  stopifnot(nrow(res) == 6L)
  stopifnot(all(c("model", "method", "delta_mu", "mean_bias", "r_kkT",
                  "feas_gap", "fidelity", "status", "f_dist",
                  "converged", "actual_miss") %in% names(res)))
  stopifnot(all(res$model == model))
  stopifnot(setequal(res$method,
                     c("FIML", "FIML_lavPredict", "MICE",
                       "missForest", "PRISM", "PRISM_MI")))
  invisible(res)
}

params_gcm <- data.frame(model = "GCM", n = 150, miss = 0.15,
                         dist = "Normal", mech = "MAR", stringsAsFactors = FALSE)
set.seed(1)
res_gcm <- run_iteration(1, params_gcm)
res_gcm <- check_schema(res_gcm, "GCM", "GCM")

# GCM MAR realized missingness follows the calibrated linear-accumulation
# formula: per-step budget m chosen so realized cellwise rate = nominal rate
m_step <- round(params_gcm$n * params_gcm$miss * t_points /
                sum((1:(t_points - 1)) * ((t_points - 1):1)))
expected_miss <- sum(vapply(1:(t_points - 1), function(t) t * m_step * (t_points - t),
                           numeric(1))) / (params_gcm$n * t_points)
cat(sprintf("GCM MAR actual_miss = %.3f (expected %.3f)\n",
            res_gcm$actual_miss[1], expected_miss))
stopifnot(abs(res_gcm$actual_miss[1] - expected_miss) < 0.01)

# regression guard: the miss=0.30 MAR condition must not degenerate into a
# fully missing final wave (the calibration keeps T4 partially observed)
params_hard <- data.frame(model = "GCM", n = 100, miss = 0.30,
                          dist = "Normal", mech = "MAR", stringsAsFactors = FALSE)
set.seed(4)
res_hard <- run_iteration(1, params_hard)
check_schema(res_hard, "GCM", "GCM hard condition (n=100, miss=0.30 MAR)")
stopifnot(res_hard$actual_miss[1] > 0.25 && res_hard$actual_miss[1] < 0.35)

# mean-targeting invariant: the in-engine joint (mu, Sigma) projection must
# pin the completed column means to the FIML model-implied means exactly
# (target_means = TRUE is the new default path); delta_mu vs the missForest
# initialiser is nonzero by design and is no longer asserted here.
prism_row <- res_gcm[res_gcm$method == "PRISM", ]
stopifnot(!is.na(prism_row$delta_mu))
stopifnot(!is.na(prism_row$mean_gap), prism_row$mean_gap < 1e-8)
stopifnot(!is.na(prism_row$lambda_sigma), prism_row$lambda_sigma == 10)
cat(sprintf("GCM PRISM mean_gap = %.2e (target_means = TRUE), delta_mu = %.2e\n",
            prism_row$mean_gap, prism_row$delta_mu))

# FIML sanity: slope variance should be near the true value 1
fiml_row <- res_gcm[res_gcm$method == "FIML", ]
cat(sprintf("GCM FIML s_var = %.3f (true 1.0)\n", fiml_row$s_var))

params_sem <- data.frame(model = "SEM", n = 150, miss = 0.15,
                         dist = "Normal", mech = "MAR", stringsAsFactors = FALSE)
set.seed(2)
res_sem <- run_iteration(1, params_sem)
res_sem <- check_schema(res_sem, "SEM", "SEM")

prism_row2 <- res_sem[res_sem$method == "PRISM", ]
stopifnot(!is.na(prism_row2$delta_mu))
stopifnot(!is.na(prism_row2$mean_gap), prism_row2$mean_gap < 1e-8)
cat(sprintf("SEM PRISM mean_gap = %.2e (target_means = TRUE), delta_mu = %.2e\n",
            prism_row2$mean_gap, prism_row2$delta_mu))

fiml_row2 <- res_sem[res_sem$method == "FIML", ]
cat(sprintf("SEM FIML b21 = %.3f (true 0.5), psi_F2 bias = %.1f%% (expect ~0)\n",
            fiml_row2$est_b21, fiml_row2$bias_var_F2))

# two-level PRISM_MI path (bootstrap-weighted forest initializer)
two_level_mi <- TRUE
set.seed(3)
res_2l <- run_iteration(1, params_gcm)
res_2l <- check_schema(res_2l, "GCM", "GCM two-level")
stopifnot(!is.na(res_2l[res_2l$method == "PRISM_MI", "r_kkT"]))
cat("GCM two-level PRISM_MI r_kkT =",
    format(res_2l[res_2l$method == "PRISM_MI", "r_kkT"], digits = 3), "\n")

cat("\nSmoke test passed.\n")
