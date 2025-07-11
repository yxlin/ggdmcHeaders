# ggdmcHeaders

The package is a collection of the 'C++' implementation of the choice
response time model. It connects the model to the Differential
Evolution Markov Chain Monte Carlo (DE-MCMC) sampler implemented in
the _ggdmc_ package.

The package supports the hierarchical modelling, Bayesian inference,
choice response time models and factorial designs, allowing users to
build their own design-based models.

The package serves as the C++ backends for the following packages:
_ggdmcModel_, _ggdmcPrior_, _ggdmcLikelihood_, _lbaModel_, _ddModel_ and _ggdmc_.

# Prerequisites
R (>= 3.5.0), Rcpp (>= 1.0.7), and RcppArmadillo 

# Installation

You can download the package from [CRAN: ggdmcHeaders](https://cran.r-project.org/web/packages/ggdmcHeaders/index.html) and install it using the source tarball.  

```
install.packages("ggdmcHeaders_0.2.9.1.tar.gz", repos = NULL, type = "source")

```

Alternatively, you can install it directly in R by running the following command:

```
install.packages("ggdmcHeaders")
```


