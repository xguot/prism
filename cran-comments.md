## Resubmission / Update (v0.1.2)
This is a patch release addressing a test failure on Windows and expanding the unit testing suite.

### Fixes
* **Windows Support:** Fixed an ERROR in tests where `mclapply` was called with `mc.cores > 1` on Windows. Added a platform check to fall back to a single core.
* **Refinement Wrappers:** Added a `robust` parameter to `smriti_forest`, `smriti_mice`, and `smriti_ranger` to allow user-controlled robust covariance projection.
* **Testing:** Expanded `tests/test_imputation.R` to include professional-grade coverage for MI, wrappers, and internal mathematical projections.
* **Simulation:** Relocated the full-grid simulation study to `hpc/` and added support for SLURM array jobs for scientific validation.

## Test Environments
* local x86_64-pc-linux-gnu (Fedora 44), R 4.6.0
* win-builder (devel)

## R CMD check results (v0.1.2)
0 errors | 0 warnings | 2 notes

*   This is a new submission.
*   Compilation used some non-portable flags (system-specific).
