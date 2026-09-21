# Shared loader for the HPC results (used by analysis.R and plots.R).
#
# Reads every per-task file sim_results/prod_results_<id>.rds, drops stale
# files written by the v1 run (old schema: 29 columns, no `model` column),
# binds the rows, and reports which task ids are missing or stale so an
# incomplete pull cannot silently enter the analysis.

suppressPackageStartupMessages({
  library(dplyr)
})

load_results <- function(dir = "sim_results", expected_ids = 1:384, warn = TRUE) {
  files <- list.files(dir, pattern = "prod_results_[0-9]+\\.rds", full.names = TRUE)
  ids <- as.integer(sub(".*prod_results_([0-9]+)\\.rds", "\\1", files))

  if (warn) {
    missing <- setdiff(expected_ids, ids)
    if (length(missing) > 0) {
      warning("Missing task files: ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  keep <- logical(length(files))
  stale <- character(0)
  for (k in seq_along(files)) {
    d <- readRDS(files[k])
    if (is.data.frame(d) && "model" %in% names(d) && nrow(d) > 0) {
      keep[k] <- TRUE
    } else {
      stale <- c(stale, basename(files[k]))
    }
  }

  if (warn && length(stale) > 0) {
    warning("Dropping ", length(stale),
            " stale v1-format file(s) (old 29-column schema): ",
            paste(stale, collapse = ", "), call. = FALSE)
  }

  res <- bind_rows(lapply(files[keep], readRDS))

  if (warn) {
    cond_gcm <- nrow(distinct(res[res$model == "GCM", c("N", "miss", "dist", "mech")]))
    cond_sem <- nrow(distinct(res[res$model == "SEM", c("N", "miss", "dist", "mech")]))
    cat(sprintf("Loaded %d new-format task files: %d rows\n", sum(keep), nrow(res)))
    cat(sprintf("Conditions represented: GCM=%d/192, SEM=%d/192\n", cond_gcm, cond_sem))
  }

  res
}
