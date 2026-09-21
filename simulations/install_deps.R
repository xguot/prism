# Install the simulation dependencies into the local Rivanna R library.
#
# Runs in a clean environment (no R_LIBS_USER set; UVA's R profile rewrites
# library paths when R_LIBS_USER is present and can break R startup).  The
# library is passed explicitly to install.packages and prepended via
# .libPaths so requireNamespace probes the right tree.

lib_path <- path.expand("~/R/rivanna-lib")
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

dependencies <- c(
  "missForest", "mice", "missRanger", "ranger",
  "lavaan", "Rcpp", "RcppArmadillo", "MASS"
)

for (pkg in dependencies) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg))
    install.packages(pkg, lib = lib_path, repos = "https://cloud.r-project.org")
  }
}

message("Dependency check complete: ", lib_path)
