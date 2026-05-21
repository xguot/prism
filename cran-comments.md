## Resubmission / Update (v0.1.1)
This is a patch release fixing a critical stability bug and aligning documentation with the implementation.

### Fixes
* **Stability:** Fixed a crash (eigenvalue error) when a column is 100% missing. Added a guard clause that stops with an informative message.
* **Consistency:** Threaded the `tol` parameter from the R interface into the C++ optimization loop (previously hardcoded to 1e-6).
* **Documentation:** Corrected `README.md` which incorrectly referenced an "MCD" estimator (the implementation uses Spearman + MAD).
* **Defaults:** Adjusted the `lambda` default to 1.0 in code and documentation to match the vignette and improve convergence out-of-the-box.
* **Testing:** Wrapped simulation tests in `requireNamespace` checks to ensure clean checks on environments without suggested packages.

## Test Environments
* local x86_64-pc-linux-gnu (Fedora 44), R 4.6.0
* Rivanna HPC cluster (CentOS Linux 7), R 4.3.1
* win-builder (devel)

## R CMD check results (v0.1.1)
0 errors | 0 warnings | 2 notes

*   This is a new submission.
*   Compilation used some non-portable flags (system-specific).
