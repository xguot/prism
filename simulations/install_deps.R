# Install dependencies for prism
lib_path <- "~/R/rivanna-lib"
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

dependencies <- c("missForest", "mice", "ranger",
                  "lavaan", "Rcpp", "RcppArmadillo", "MASS")

for (pkg in dependencies) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg))
    install.packages(pkg, lib = lib_path, repos = "https://cloud.r-project.org")
  }
}
