library(testthat)
library(prism)

# R CMD check runs this file with tests/ as the working directory, where
# test_check() resolves the testthat/ tree.  Standalone runs from the
# package root (R_LIBS=libs Rscript tests/testthat.R) need the explicit
# path instead.
if (dir.exists("testthat")) {
  test_check("prism")
} else {
  # prism is already attached via library() above
  test_dir("tests/testthat", load_package = "none")
}
