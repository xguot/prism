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

# mean-preservation invariant: completed PRISM data must not drift from X0
prism_row <- res_gcm[res_gcm$method == "PRISM", ]
stopifnot(!is.na(prism_row$delta_mu))
cat(sprintf("GCM PRISM delta_mu = %.2e (expect ~0)\n", prism_row$delta_mu))
stopifnot(prism_row$delta_mu < 1e-8)

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
cat(sprintf("SEM PRISM delta_mu = %.2e (expect ~0)\n", prism_row2$delta_mu))
stopifnot(prism_row2$delta_mu < 1e-8)

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
