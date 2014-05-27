#ifndef TestProblems_h
#define TestProblems_h

class TH1D;
class TH2D;
class TF1;

#include <TMatrixD.h>
#include <TVectorD.h>

struct TestProblem       // Output from test problem generator.
{
  TMatrixD R;            // Fine-binned response matrix.
  TH2D *Response;        // Response matrix.
  TH1D *xTruth;          // Discrete version of true PDF.
  TH1D *xTruthEst;       // Estimator for xTruth.
  TH1D *xIni;            // Model truth.
  TH1D *bIdeal;          // Observed b, no noise.
  TH1D *bNoisy;          // Observed b, perturbed by noise.
  TH1D *eff;             // Efficiency: meas/true vs. true x
};

// Test problems
// Shaw: sets up discretized n x n kernel A, as
// well as true and measured vectors x and b (both size n).
// C. B. Shaw, Jr., "Improvements of the resolution of an instrument
// by numerical solution of an integral equation", J. Math. Anal. Appl. 37
// (1972), 83–112.
TestProblem ShawSystem(const int n, double noise=0.);
void ShawSystem(const int n, TMatrixD &A, TVectorD &x, TVectorD &b,
                double noise=0);

TestProblem MonteCarloConvolution(const int m,
                                  const int n,
                                  const double xm1,
                                  const double xm2,
                                  const double xt1,
                                  const double xt2,
                                  TF1 *truthFn,
                                  TF1 *kernelFn,
                                  const int nEvents);

TestProblem AtlasDiJetMass(const int Nt,
                           const int Nr,
                           double tbins[],
                           double rbins[],
                           const double apar = 0.5,
                           const double bpar = 0.1,
                           const int nEvents = int(1e7),
                           const double evtWeight = 1e-3);
#endif
