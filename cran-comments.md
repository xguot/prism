## Resubmission Notes
This is a resubmission. The following algorithmic issues from internal review have been addressed:

### 1. Gradient Descent (src/lagrange.cpp)
- **Corrected the gradient derivation.** The covariance-projection gradient now uses column-centered data (`X_tilde`) rather than raw `X`, preventing unintended location shifts from leaking into the covariance constraint.
- **Added fidelity term.** The loss function is now `L(X) = (1/2)||X - X_imp||² + (lambda/2)||cov(X) - Sigma_target||²`, anchoring the solution near the initial imputation. The previous code minimised only `||cov(X) - Sigma_target||²`, allowing the imputed values to drift arbitrarily far from the original.
- **Numerical safety.** Added `has_nan()` / `has_inf()` guards that abort with an informative error if optimisation diverges (typically due to an ill-conditioned target matrix).
- **Convergence check.** Moved the Frobenius-norm tolerance check to the *post-update* state so that convergence is detected on the same iteration it occurs, not one iteration late.

### 2. Robust Target Covariance (R/router.R)
- **PSD projection.** Pairwise-deletion correlation/covariance matrices are not guaranteed positive semidefinite. A `nearest_psd()` helper (Higham 1988) now projects both the Spearman correlation matrix and the final `D*R*D` matrix onto the PSD cone before they are passed to the C++ backend.
- **Eigenvalue validation.** Added an explicit `eigen()` check that aborts with an informative error if the target matrix is numerically singular, instead of silently feeding a degenerate target to the optimiser.
- **Documentation.** The robust path is now accurately described as "pairwise Spearman + MAD → nearest PSD" rather than the previously claimed (but unimplemented) MCD estimator.

### 3. Multiple Imputation (R/smriti_mi.R)
- **PSD safeguard.** In addition to `anyNA()`, the bootstrap-sampled correlation matrix is now checked for positive semidefiniteness before the sample is accepted.
- **Error resilience.** The `smriti_impute()` call is wrapped in `tryCatch()` so that a singular-target failure in one bootstrap draw retries rather than aborting the entire MI loop.
- **Documentation.** Added a prominent `@details` section explaining that this is *approximate* bootstrap MI (not proper Bayesian MI) and that standard errors may be mildly anti-conservative.

### 4. Hyperparameter Defaults
- `lambda`: 0.5 → **1.0** (balanced fidelity-constraint trade-off with both loss terms)
- `learning_rate`: 1e-7 → **0.001** (effective step size with the augmented loss)
- `max_iter`: 1000 → **2000** (sufficient for convergence with practical step sizes)
- Added a post-optimisation convergence diagnostic that warns if the final Frobenius distance exceeds `tol * 100`.

## Test Environments
* local x86_64-pc-linux-gnu (AMD Ryzen), R 4.3.1
* Rivanna HPC cluster (CentOS Linux 7), R 4.3.1
* win-builder (devel)

## R CMD check results
0 errors | 0 warnings | 0 notes
