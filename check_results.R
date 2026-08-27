files <- list.files("sim_results", pattern = "prod_results_[0-9]+\\.rds", full.names = TRUE)
ids <- as.integer(sub(".*prod_results_([0-9]+).rds", "\\1", files))
cat("task files:", length(files), "\n")
cat("duplicate ids:", paste(ids[duplicated(ids)], collapse = ","), "\n")
cat("missing ids:", if (length(setdiff(1:384, ids))) paste(setdiff(1:384, ids), collapse = ",") else "NONE", "\n")
for (f in files) {
  d <- readRDS(f)
  ok <- is.data.frame(d) && "model" %in% names(d) && nrow(d) > 0
  cat(sprintf("%-28s rows=%-5d cols=%-3d new_schema=%s\n",
              basename(f), if (is.data.frame(d)) nrow(d) else 0L,
              if (is.data.frame(d)) ncol(d) else 0L, ok))
}
