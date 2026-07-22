# lagrange

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**lagrange** is an R package that produces a completed longitudinal dataset whose covariance structure matches the model-implied covariance from a FIML-estimated latent growth model. It fits the model via `lavaan`, extracts the structural target, and projects an initial imputation onto that manifold using a C++ Lagrangian-constrained gradient descent engine.

## Installation

```R
# Development version
# install.packages("devtools")
devtools::install_github("xguot/lagrange")
```

## Usage

```R
library(lagrange)
library(lavaan)

# Define a latent growth model
model <- "
  i =~ 1*T1 + 1*T2 + 1*T3 + 1*T4
  s =~ 0*T1 + 1*T2 + 2*T3 + 3*T4
"

# Impute missing values while preserving the FIML-implied structure
imputed_data <- lagrange_fiml(
  data   = clinical_df,
  model  = model,
  lambda = 0.5
)
```

## Architecture

The pipeline executes in two phases:

1. **FIML Estimation:** A latent growth model is fit to the incomplete data via full-information maximum likelihood (`lavaan::growth`). The model-implied covariance matrix Σ_FIML is extracted as the structural target.
2. **Lagrangian Projection:** An initial imputation (user-supplied or column-mean fallback) is projected onto the Σ_FIML manifold via constrained gradient descent. Only originally-missing cells are updated; observed data is held fixed.

The C++ engine minimizes:

$$L(X) = \|\operatorname{cov}(X) - \Sigma_{\text{FIML}}\|_F^2$$

subject to observed cells frozen in place.

## Why lagrange

FIML estimates parameters correctly under MAR but does not fill in missing values. `lavPredict(type="ov")` (conditional expectations) attenuates variance by 40–60%. lagrange produces a completed dataset whose structural parameters match the FIML model without variance attenuation — enabling downstream analyses that require complete data.

## Citation

If you use **lagrange** in your research, please cite:

> Guo, X. (2026). lagrange: FIML Covariance Projection for Longitudinal Missing Data. R package version 0.2.0.
