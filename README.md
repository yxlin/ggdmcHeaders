# ggdmcHeaders

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ggdmcHeaders)](https://cran.r-project.org/package=ggdmcHeaders)
[![Downloads](https://cranlogs.r-pkg.org/badges/ggdmcHeaders)](https://cran.r-project.org/package=ggdmcHeaders)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/yxlin/ggdmcHeaders/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yxlin/ggdmcHeaders/actions/workflows/R-CMD-check.yaml)


**ggdmcHeaders 0.2.9.2 (development)** provides the C++ backend for choice response time models, enabling high-performance computation for hierarchical modelling and Bayesian inference.  

This version introduces **Cognitive Diagnostic Model (CDM)** support, extending the package beyond LBA and DDM to provide a unified backend for a wider class of cognitive models.

It integrates with the Differential Evolution Markov Chain Monte Carlo (DE-MCMC) sampler from the [`ggdmc`](https://cran.r-project.org/package=ggdmc) package, allowing users to build and fit design-based cognitive models, including factorial designs.

This package serves as a shared C++ codebase for the following packages:  

- [`ggdmcModel`](https://cran.r-project.org/package=ggdmcModel)
- [`ggdmcPrior`](https://cran.r-project.org/package=ggdmcPrior)
- [`ggdmcLikelihood`](https://cran.r-project.org/package=ggdmcLikelihood)
- [`lbaModel`](https://cran.r-project.org/package=lbaModel)
- [`ddModel`](https://cran.r-project.org/package=ddModel)
- [`ggdmc`](https://cran.r-project.org/package=ggdmc)
- `cdModel` – Cognitive Diagnostic Models (new integration)

---

## ✨ Features
- C++ implementations of choice response time models.
- Full support for hierarchical modelling and Bayesian inference.
- Optimised for DE-MCMC sampling in `ggdmc`.
- Works seamlessly with multiple model packages in the `ggdmc` ecosystem.
- Designed for flexibility in building **design-based cognitive models**.
- **New in 0.2.9.2**: Support for **Cognitive Diagnostic Models (CDM)**, enabling latent skill diagnosis and mastery modelling.

---

## 📦 Prerequisites
- **R** ≥ 3.5.0  
- **Rcpp** ≥ 1.0.7  
- **RcppArmadillo**  

---

## 📥 Installation

### Development version from GitHub

```r
# install.packages("remotes")
remotes::install_github("yxlin/ggdmcHeaders", ref = "dev")
```


## 🔗 Related Packages
- ggdmc – Main modelling framework
- ggdmcModel – Model specification
- ggdmcPrior – Prior distributions
- ggdmcLikelihood – Likelihood evaluation
- lbaModel – Linear Ballistic Accumulator models
- ddModel – Diffusion Decision models
- cdModel – Cognitive Diagnostic Models (new integration)

## 📄 License
This package is released under the GPL (≥ 3) license.