# smriti: Structural Variance Preservation via Lagrangian Manifold Routing

[![CRAN status](https://www.r-pkg.org/badges/version/smriti)](https://CRAN.R-project.org/package=smriti)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**smriti** is an R package for automated longitudinal missing data imputation. Unlike standard non-parametric methods that may drift from the underlying data structure, **smriti** utilizes a Lagrangian constraint mechanism to project machine-learning hallucinations back toward a target covariance manifold.

## Key Features

- **Structural Integrity:** Ensures that imputed values respect the established covariance structure of the observed data.
- **Robust Estimation:** Optional MCD (Minimum Covariance Determinant) estimator to suppress the influence of outliers on the target manifold.
- **High Performance:** Core Lagrangian projection engine implemented in C++ via Rcpp and Armadillo for maximum efficiency on large longitudinal datasets.
- **Hybrid Approach:** Combines the flexibility of Random Forest-based imputation (`missForest`) with the rigor of structural manifold routing.

## Installation

You can install the released version of **smriti** from [CRAN](https://CRAN.R-project.org) with:

```R
install.packages("smriti")
```

Or the development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("xguot/smriti")
```

## Quick Start

```R
library(smriti)

# Example: Impute missing values in a longitudinal dataset
# data: Your data frame
# time_cols: The columns representing longitudinal measurements
imputed_data <- smriti_impute(
  data = clinical_df, 
  time_cols = c("V1", "V2", "V3", "V4"),
  lambda = 0.5,
  robust = TRUE
)
```

## Methodology

The package operates in three phases:
1. **Initial Hallucination:** Uses a non-parametric approach to generate initial dense point-clouds.
2. **Manifold Mapping:** Establishes a target covariance manifold from observed data (optionally using robust estimators).
3. **Lagrangian Routing:** Projects trajectories back to the structural manifold using a gradient descent update constrained by the penalty weight $\lambda$.

## Citation

If you use **smriti** in your research, please cite:

> Guo, X. (2026). smriti: Structural Variance Preservation for Longitudinal Missing Data Imputation. R package version 0.1.0.

*Manuscript in preparation.*

## License

This package is licensed under the MIT License.
