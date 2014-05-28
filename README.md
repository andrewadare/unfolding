# unfolding.git

### C++/ROOT package for unfolding and inference

This package provides many tools for tackling discrete linear inverse problems, including

- Singular value decomposition and Generalized SVD
    + Extension of the method described in [SVD approach to data unfolding](http://arxiv.org/abs/hep-ph/9509307)
    + GSVD provides more regularization flexibility than the classic SVD method, and allows rectangular (m > n) matrices.
    + Error propagation through full covariance matrix
- Preconditioned conjugate gradients for least squares (P-CGLS)
    + Fast iterative method for large problems allowing flexible regularization - see [Discrete Inverse Problems](http://epubs.siam.org/doi/book/10.1137/1.9780898718836) by P. C. Hansen
- Richardson-Lucy iterative solver
    + [Bayesian-based iterative method of image restoration](http://www.opticsinfobase.org/josa/abstract.cfm?id=54565) by W. H. Richardson
    + See also "An iterative technique for the rectification of observed distributions" L. B. Lucy, The Astronomical Journal, vol 79, no 6, 1974.
- $\chi^2$ minimizer using ROOT's ```TMinuit2``` package (requires ROOT to be built with ```--enable-minuit2```)
- "Fully Bayesian" method using MCMC Metropolis random walk sampler.
    + ["Fully Bayesian Unfolding"](http://arxiv.org/pdf/1201.4612v4.pdf) by  G. Choudalakis

Significant attention has been devoted to visualization. The examples produce spectral decomposition plots, step-wise animations for the iterative solvers, and surface plots demonstrating evolution the solutions with regularization strength.

Several quantitative methods are included to assist in finding the optimal regularization strength, including generalized cross-validation, L-curve analysis, effective degrees of freedom, and singular value spectra.

### Getting started
There is no build system other than ROOT's ACLiC compiler (see ```examples/rootlogon.C```). This package requires ```MatrixUtils.h``` and ```UtilFns.h``` from [andrewadare/utils](https://github.com/andrewadare/utils). Either symlink them to the same directory as ```README.md```, or add their location to the ACLiC and interpreter include paths by editing ```examples/rootlogon.C```. The example scripts should then work without further setup; for example do ```$ root ConvolutionExample.C+```. 

This code has been tested on Mac OS X 10.9 against ROOT 5 (5.34/18) and ROOT 6 beta (5.99/06). Unfortunately, I can't say much about any other setup.
