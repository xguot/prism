# prism

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**prism** is an R package that produces a completed dataset whose covariance structure matches the model-implied covariance from a FIML-estimated structural equation model (latent growth models, CFAs, and general SEMs). It fits the model via `lavaan`, extracts the structural target, and projects an initial imputation onto that target under a regularized dual objective: fidelity to the initial imputation (nonparametric, local information) is traded off against covariance matching (parametric structure) through `lambda_sigma`, while the column means of the initial imputation are preserved exactly. A C++ engine performs mean-preserving projected gradient descent with KKT convergence certification; `prism_mi` adds two-level proper multiple imputation.

## Installation

```R
# Development version
# install.packages("devtools")
devtools::install_github("xguot/prism")
```

## Usage

```R
library(prism)
library(lavaan)

# Define a latent growth model
model <- "
  i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
  s =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
"

# Impute missing values while preserving the FIML-implied structure
imputed_data <- prism_fiml(
  data         = clinical_df,
  model        = model,
  lambda_sigma = 0.5
)

# Convergence and constraint diagnostics
attr(imputed_data, "prism_diagnostics")

# General SEM support: any single-group lavaan model (CFA, structural
# regressions, mediation) runs through the same engine; model syntax or a
# pre-fitted lavaan object are both accepted
sem_result <- prism_sem(
  data  = sem_df,
  model = "eta1 =~ y1 + y2 + y3
           eta2 =~ y4 + y5 + y6
           eta2 ~ eta1"
)

# Proper multiple imputation (two-level: stochastic initializations +
# parameter draws)
mi_data <- prism_mi(
  data    = clinical_df,
  fit     = fit,
  m       = 20,
  seed    = 42,
  lambda_sigma = 0.5
)
```

## Architecture

The pipeline executes in two phases:

1. **FIML Estimation:** A structural equation model (growth model, CFA, or general SEM) is fit to the incomplete data via full-information maximum likelihood (`lavaan::growth` or `lavaan::sem`). The model-implied covariance matrix Σ_FIML of the observed variables is extracted as the structural target.
2. **Regularized Projection:** An initial imputation (user-supplied or column-mean fallback) is projected toward the Σ_FIML target under the dual objective

$$f(X) = \frac{1}{2}\|M \odot (X - X^{(0)})\|_F^2 + \frac{\lambda_\Sigma}{2}\|\operatorname{Cov}(X) - \Sigma_{\text{FIML}}\|_F^2$$

subject to the observed cells being frozen and the column means of X⁽⁰⁾ preserved. `lambda_sigma` is a genuine structural trade-off (it changes the fixed points, not the step size): `0` returns the initial imputation, large values enforce the FIML structure. The C++ engine computes the analytic gradient, centers the missing-cell gradient at zero within each column before every step (mean-preserving projected gradient descent), and uses Armijo backtracking on the full objective.

**Convergence certification.** The engine reports KKT diagnostics rather than only a covariance distance: the projected-gradient norm `r_kkT` (stationarity), the covariance feasibility gap `feas_gap`, and per-column Lagrange multiplier estimates `nu`. The `status` field distinguishes a feasible optimum (`converged_feasible`) from a constrained geometric limit with binding mean constraints (`constrained_geometric_limit`), an unreachable target (`geometric_infeasible`), and non-stationary stops (`stalled_line_search`, `max_iter_reached`).

## Scope

`prism_sem()` supports any single-group, single-level SEM with continuous observed variables. Ordinal (WLSMV), multi-group, and multilevel models are rejected with clear errors; multi-group data can be projected group by group. Model-implied means are deliberately not targeted: first moments are anchored to the initial imputation.

## Multiple Imputation

`prism_mi` implements two-level proper multiple imputation. Level 1 injects missing-value (predictive) uncertainty by drawing a fresh stochastic initial imputation for each completed dataset — by default a bootstrap-weighted `missRanger` forest; any stochastic imputer can be supplied via `initializer`. Level 2 draws perturbed model parameters θ⁽ᵏ⁾ ~ N(θ̂, ACOV(θ̂)) and reconstructs Σ⁽ᵏ⁾. Both sources of randomness enter the between-imputation variance under Rubin's rules.

## Why prism

FIML estimates parameters correctly under MAR but does not fill in missing values. `lavPredict(type="ov")` (conditional expectations) attenuates variance by 40–60%. prism produces a completed dataset whose structural parameters match the FIML model without variance attenuation — enabling downstream analyses that require complete data.

## Citation

If you use **prism** in your research, please cite:

> Guo, X. (2026). prism: FIML Covariance Projection for Longitudinal Missing Data. R package version 0.3.0.
