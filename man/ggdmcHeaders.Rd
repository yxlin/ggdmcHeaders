% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggdmcHeaders.R
\name{ggdmcHeaders}
\alias{ggdmcHeaders}
\title{'C++' Backend for the 'ggdmc' Ecosystem}
\description{
The package is a collection of the 'C++' implementation of the choice
response time model. It connects the model to the Differential
Evolution Markov Chain Monte Carlo (DE-MCMC) sampler implemented in
the \pkg{ggdmc} package.
}
\details{
The package supports the hierarchical modelling, Bayesian inference,
choice response time models and factorial designs, allowing users to
build their own design-based models.

The package serves as the C++ backends for the following packages:
\pkg{ggdmcModel}, \pkg{ggdmcPrior}, \pkg{ggdmcLikelihood},
\pkg{lbaModel}, '\pkg{ddModel} and \pkg{ggdmc}.
}
\references{
\itemize{
   \item Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting
         the parameters of the diffusion model: A empirical validation
         \emph{Memory and Cognition}, 32(7), 1206--1220.
   \item Ratcliff, R., & McKoon, G. (2008). The diffusion decision
         model: Theory and data for two-choice decision tasks.
         \emph{Neural Computation}, 20(4), 873--922.
   \item Brown S., & Heathcote, A. (2008). The simplest complete model
         of choice response time: Linear ballistic accumulation.
         \emph{Cognitive Psychology}, 57(3), 153--178.
 }
}
