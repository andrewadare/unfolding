// ObjectiveFns.h
// Stateful functions for use in optimization algorithms.

#include <TMath.h>
#include <TVectorD.h>
#include <vector>
#include <iostream>

using namespace std;

double PoissonLikelihood(const TVectorD &b, const TVectorD &R);
double LogPoissonLikelihood(const TVectorD &b, const TVectorD &R);
double LogPoisson(double x, double mu);
double LogGaussian(double x, double mu, double sigma, bool norm);
double LogFactorial(int n);

// Base -----------------------------------------------------------------------

// Log likelihood function base functor
// This is specifically intended for a linear system Ax = b, where x is a
// trial solution, b is a data vector, and A is the transfer matrix.
// The likelihood is thus computed from a comparison of Ax and b.
struct LogLikeFn
{
  explicit LogLikeFn(TMatrixD &A, TVectorD &b) : fMatrix(A), fData(b) {}
  virtual double operator()(const TVectorD & /* x */) = 0;
  virtual ~LogLikeFn() {}

  TMatrixD fMatrix;
  TVectorD fData;
};

// Log Likelihood function allowing multiple model / dataset pairs.
struct LLMultiFn
{
  explicit LLMultiFn(vector<TMatrixD> &avec, vector<TVectorD> &bvec,
                     vector<double> &weights) :
    fAVec(avec), fbVec(bvec), fwVec(weights) {}
  virtual double operator()(const TVectorD & /* x */) = 0;
  virtual ~LLMultiFn() {}

  vector<TMatrixD> fAVec;
  vector<TVectorD> fbVec;
  vector<double>   fwVec;
};

// Prior function base functor
struct LogPrior
{
  explicit LogPrior(vector<double> &pars) : fPars(pars) {}

  virtual double operator()(const TVectorD & /* trialpars */) = 0;
  virtual ~LogPrior() {}

  vector<double> fPars; // Prior hyperparameters
};

// Inherited ------------------------------------------------------------------

// Log Poisson likelihood fn. for linear system A*x = b
struct LogPoissonLikeFn : public LogLikeFn
{
  explicit LogPoissonLikeFn(TMatrixD &A, TVectorD &b) : LogLikeFn(A, b) {}

  double operator()(const TVectorD &x)
  {
    TVectorD Ax = fMatrix*x;
    return LogPoissonLikelihood(fData, Ax);
  }
};

// Log Poisson likelihood fn. for multiple linear systems A[i]*x = b[i]
struct PoissonLLMultiFn : public LLMultiFn
{
  explicit PoissonLLMultiFn(vector<TMatrixD> &avec,
                            vector<TVectorD> &bvec,
                            vector<double> &weights) :
    LLMultiFn(avec, bvec, weights) {}

  double operator()(const TVectorD &x)
  {
    double ll = 0;
    int nsets = fAVec.size();
    for (int i=0; i<nsets; i++)
    {
      TVectorD Ax = fAVec[i]*x;
      ll += fwVec[i]*LogPoissonLikelihood(fbVec[i], Ax);
    }
    return ll;
  }
};

// Log Gaussian prior. Penalize disagreement between
// solution and fPars[1..Nt]. Note offset from first entry:
// alpha = fPars[0] is the prior precision (reciprocal variance).
// Use this, for example, to bias your result towards the model truth.
// Returns a positive quantity.
/*
// Example setup:
std::vector<double> regpars(Nt+1, 0.0); // Prior precision [0], MC truth [1+]
regpars[0] = 1.0;
for (int i=0; i<Nt; i++)
  regpars[i+1] = T(i);
Chi2RegFn regfunc(regpars);
*/

struct Chi2RegFn : public LogPrior
{
  explicit Chi2RegFn(vector<double> &pars) : LogPrior(pars) {}

  double operator()(const TVectorD &trialpars)
  {
    double alpha = fPars[0], mu = 0., chi2 = 0.;
    for (int t=0; t<trialpars.GetNrows(); t++)
    {
      mu = fPars[t+1];
      // printf(" %.0f", trialpars(t));
      double arg = (trialpars(t) - mu)/mu;
      chi2 += 0.5*alpha*arg*arg;
    }
    return chi2;
  }
};

// Total curvature. Discontinuity excluded at i = iskip if requested.
// Eq. (38), NIM A 372 (1996) 469-481.
struct CurvatureRegFn : public LogPrior
{
  explicit CurvatureRegFn(vector<double> &pars) : LogPrior(pars) {}

  double operator()(const TVectorD &x)
  {
    double alpha = fPars[0];
    int iskip    = fPars.size() > 1 ? fPars[1] : -1;
    double delta = 0;
    double curv  = 0;

    for (int i=1; i<x.GetNrows()-1; i++)
    {
      if (x(i) == 0)
        continue;

      if (i == iskip - 1)
        delta = x(i-1) - x(i);
      else if (i == iskip)
        delta = x(i+1) - x(i);
      else
        delta = x(i+1) - 2*x(i) - x(i-1);

      curv += delta*delta / x(i);
    }
    return alpha*alpha*TMath::Log(curv);
  }
};

// Supporting functions -------------------------------------------------------
double
PoissonLikelihood(const TVectorD &b, const TVectorD &R)
{
  double ans = 1;
  for (int r=0; r<b.GetNrows(); ++r)
  {
    ans *= TMath::Poisson(b(r), R(r));
  }
  return ans;
}

double
LogPoissonLikelihood(const TVectorD &b, const TVectorD &R)
{
  double ans = 0;
  for (int r=0; r<b.GetNrows(); ++r)
  {
    ans +=LogPoisson(b(r), R(r));
  }
  return ans;
}

double
LogPoisson(double x, double mu)
{
  if (mu >= 1000)
    return LogGaussian(x,mu,TMath::Sqrt(mu),true);

  if (x < 0)
    return 0;

  if (x==0.)
    return -mu;

  return x*TMath::Log(mu) - mu - LogFactorial(int(x));
}

double
LogGaussian(double x, double mu, double sigma, bool norm)
{
  if (sigma <= 0)
    return 0;

  double ans = -0.5*(x-mu)*(x-mu)/sigma/sigma;

  if (norm)
    ans -= TMath::Log(TMath::Sqrt(TMath::TwoPi())*sigma);

  return ans;
}

double
LogFactorial(int n)
{
  // Cache ln(n!) values for speed.
  static int nStoredLnFac = 10000;
  static double *lnFactorial = 0;

  if (!lnFactorial)
  {
    // Compute & store ln n! on initial function call
    if (0)
    {
      cout << Form("Pre-caching ln(n!) values up to n = %d...",
                   nStoredLnFac)
           << flush;
    }

    lnFactorial = new double[nStoredLnFac+1];

    // Set ln(0) = 0 for stability
    lnFactorial[0] = 0.;
    for (int i=1; i<=nStoredLnFac; i++)
      lnFactorial[i] = TMath::LnGamma(i+1);

    //    Printf("Done.");
  }

  if (n<=nStoredLnFac)
    return lnFactorial[n];

  return  TMath::LnGamma(n+1);
}
