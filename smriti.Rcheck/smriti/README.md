# smriti

[![CRAN status](https://www.r-pkg.org/badges/version/smriti)](https://CRAN.R-project.org/package=smriti) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**smriti** is an R package for automated longitudinal missing data imputation. It combines the predictive flexibility of non-parametric machine learning with a C++ Lagrangian projection engine to strictly preserve the structural variance of the target covariance manifold.

## Installation

```R
# Stable CRAN release
install.packages("smriti")

# Development version
# install.packages("devtools")
devtools::install_github("xguot/smriti")
```

## Usage

Impute longitudinal missing data while preserving the underlying covariance structure:

```R
library(smriti)

imputed_data <- smriti_impute(
  data = clinical_df, 
  time_cols = c("V1", "V2", "V3", "V4"),
  lambda = 0.5,
  robust = TRUE  # Enables the MCD estimator to suppress outliers
)
```

## Architecture

The imputation pipeline executes in three phases:

1. **Initialization:** Generates a dense preliminary point-cloud via Random Forest (missForest).
2. **Manifold Mapping:** Establishes the target covariance structure from observed data, with optional robust estimation.
3. **Lagrangian Routing:** Projects the initial matrix back onto the structural manifold via a constrained gradient descent update.

## Citation

If you utilize **smriti** in your research, please cite:

> Guo, X. (2026). smriti: Structural Variance Preservation for Longitudinal Missing Data Imputation. R package version 0.1.0.
