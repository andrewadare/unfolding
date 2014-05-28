# Unfolding

## Libraries and example macros to solve unfolding problems using ROOT

This package provides many tools for tackling discrete linear inverse problems, including

- Singular value decomposition and Generalized SVD
- Preconditioned conjugate gradients for least squares (P-CGLS)
- Richardson-Lucy iterative solver
- $\chi^2$ minimizer using TMinuit2 package (requires ROOT to be built with ```--enable-minuit2```)
- "Fully Bayesian" method using Markov Chain Monte Carlo. (Metropolis random walk sampler).

Significant attention has been devoted to visualization. The examples produce spectral decomposition plots, step-wise animations for the iterative solvers, and surface plots demonstrating evolution the solutions with regularization strength.

Several quantitative methods are included to assist in finding the optimal regularization strength, including generalized cross-validation, L-curve analysis, effective degrees of freedom, and singular value spectra.

This code has been tested on Mac OS X 10.9 against ROOT 5 (5.34/18) and ROOT 6 beta (5.99/06). Unfortunately, I can't say much about any other setup.
