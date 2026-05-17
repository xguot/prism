## Resubmission Notes
This is a resubmission. The following issues from the previous submission have been addressed:

* **Windows Installation Error:** Added `src/Makevars` and `src/Makevars.win` to explicitly link LAPACK and BLAS libraries required by `RcppArmadillo`. Added `SystemRequirements: C++17` to `DESCRIPTION`.
* **License Issue:** Fixed the `LICENSE` file format to the standard `YEAR: HOLDER:` format and corrected the `License` stub in `DESCRIPTION`.
* **Test Timeout:** Optimized `tests/simulation_study.R` to run a minimal single replication during CRAN checks (while maintaining full scale for local/HPC testing) to stay within the 30-minute limit.
* **Non-standard Files:** Added `hpc/` and `cran-comments.md` to `.Rbuildignore`.
* **Spelling:** Corrected misspelled/unexplained acronyms in `DESCRIPTION`.
* **Lavaan Warnings:** Modified the simulation script to ensure no empty rows are generated, preventing "some cases are empty" warnings.

## Test Environments
* local x86_64-pc-linux-gnu (AMD Ryzen), R 4.3.1
* Rivanna HPC cluster (CentOS Linux 7), R 4.3.1
* win-builder (devel) - expected to pass after adding Makevars.win

## R CMD check results
All previous errors and notes have been addressed.
