// TestProblems.C
// Self-contained examples to test and study unfolding methods.
//
// - Shaw
// - MC Convolution
// - ATLAS Dijet Mass
//
// Please report issues to Andrew Adare andrewadare@gmail.com

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TString.h>
#include <iostream>

#include "MatrixUtils.C"

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

TestProblem
MonteCarloConvolution(const int m,
                      const int n,
                      const double xm1,
                      const double xm2,
                      const double xt1,
                      const double xt2,
                      TF1 *truthFn,
                      TF1 *kernelFn,
                      const int nEvents)
{
  static int id = 0; id++;
  TestProblem t;
  double dm = (xm2-xm1)/m;
  double dt = (xt2-xt1)/n;
  TMatrixD R(m,n);

  // Discretize the kernel to fill R
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      R(i,j) = kernelFn->Eval(xt1+j*dt-i*dm);
  t.Response = MatrixUtils::Matrix2Hist(R, Form("R%d",id), xm1, xm2, xt1, xt2);
  t.Response->SetTitle(Form("%d x %d convolution matrix;"
                            "s (observed);t (true)", m, n));

  TH2D *RMC = new TH2D(Form("R_MC%d",id), "A_{MC}",
                       m,xm1,xm2,n,xt1,xt2);

  // Model a true and a measured distribution
  // There is no bIdeal for this problem.
  t.xTruth    = new TH1D("hTrue",    "",         n, xt1, xt2);
  t.xTruthEst = new TH1D("hTrueEst", "hTrueEst", n, xt1, xt2);
  t.bNoisy    = new TH1D("hMeas",    "hMeas",    m, xm1, xm2);
  t.xTruthEst->Sumw2();
  t.bNoisy->Sumw2();


  for (Int_t i=0; i<nEvents; i++)
  {
    Double_t xt = truthFn->GetRandom();
    t.xTruthEst->Fill(xt);
    Double_t xm = xt + kernelFn->GetRandom();
    t.bNoisy->Fill(xm);
    RMC->Fill(xm,xt);
  }

  // MC efficiency histogram
  t.eff = RMC->ProjectionY("eff",1,m);
  t.eff->Divide(t.xTruthEst);

  // Exact truth histogram
  for (int j=1; j<=n; j++)
  {
    double val = truthFn->Eval(t.xTruth->GetBinCenter(j));
    t.xTruth->SetBinContent(j, val);
  }
  t.xTruth->Scale(nEvents/t.xTruth->Integral());

  return t;
}

TestProblem
AtlasDiJetMass(const int Nt,
               const int Nr,
               double tbins[],
               double rbins[],
               const double apar,
               const double bpar,
               const int nEvents,
               const double evtWeight)
{
  // From "Fully Bayesian Unfolding" by G. Choudalakis.
  // see arXiv:1201.4612v4
  // Recommended binning:
  // double bins[Nt+1] = {0};
  // for (int j=0; j<=Nt; j++)
  //   bins[j] = 500*TMath::Exp(0.15*j);

  static int id = 0; id++;
  TestProblem t;
  TRandom3 ran;

  // Mass distribution dN/dM
  TF1 *mass = new TF1("mass_dist",
                      "TMath::Power(1.-x/7000,6.0)/TMath::Power(x/7000,4.8)",
                      tbins[0], tbins[Nt]);

  // Generate test problem by MC convolution
  TH1D *hT   = new TH1D("hT",   "Truth mass dist. #hat{T}", Nt, tbins);
  TH1D *hTmc = new TH1D("hTmc", "MC Truth mass dist. #tilde{T}", Nt, tbins);
  TH1D *hD   = new TH1D("hD", "Measured mass dist.", Nr, rbins);
  TH2D *hM   = new TH2D("hM", "Migration matrix", Nr, rbins, Nt, tbins);
  hM->SetTitle(Form("%d x %d migration matrix M_{tr};"
                    "mass (observed);mass (true)", Nr, Nt));
  hT->SetLineWidth(2);
  hD->SetLineWidth(2);

  std::cout << Form("Generating test problem...") << std::flush;
  // Fill histos with MC events.
  // The response matrix gets more statistics than the data.
  double xt, xm, sigma;
  for (int i=0; i<nEvents; i++)
  {
    xt = mass->GetRandom();
    sigma = apar*TMath::Sqrt(xt) + bpar*xt;
    xm = xt + ran.Gaus(0, sigma);
    hM->Fill(xm, xt);
    hTmc->Fill(xt, evtWeight);
    hD->Fill(xm, evtWeight);
  }

  // Simulate Poisson fluctuations in real data (integer content, empty bins)
  for (int r=1; r<=Nr; r++)
    hD->SetBinContent(r, ran.Poisson(hD->GetBinContent(r)));

  // The true truth \hat{T}
  double totmass = mass->Integral(tbins[0],tbins[Nt]);
  for (int j=1; j<=Nt; j++)
  {
    hT->SetBinContent(j, evtWeight*nEvents*mass->Integral(tbins[j-1],
                      tbins[j])/totmass);
    hT->SetBinError(j, 0);
  }
  cout << "Done." << endl;

  // Projection of migration matrix to truth axis. hMt is the
  // numerator for efficiency. Bin contents should be counts
  // here. Elements are normalized to contain probabilities only after
  // projection & division by T-tilde.
  TH1D *hMt  = hM->ProjectionY("hMt",1,Nr);
  TH1D *heff = (TH1D *)hMt->Clone("heff");
  heff->Divide(hTmc);
  heff->Scale(evtWeight);
  hM->Scale(1./nEvents);
  hMt->Scale(1./nEvents);

  t.Response  = hM;
  t.xTruth    = hT;
  t.xTruthEst = hTmc;
  t.xIni      = hMt;
  t.bNoisy    = hD;
  t.eff       = heff;

  return t;
}

TestProblem
ShawSystem(const int n, double noise)
{
  TMatrixD A(n,n);
  TVectorD x(n);
  TVectorD b_ideal(n);
  ShawSystem(n, A, x, b_ideal, 0.);
  TVectorD b(b_ideal);

  // Add Gaussian white noise to b
  if (noise > 0.)
  {
    TRandom3 r3;
    for (int j=0; j<n; j++) b(j) += noise*r3.Gaus();
  }
  // Assign output struct members.
  // No xTruthEst histogram for this problem. Don't ask for it!
  TestProblem t;
  t.Response = MatrixUtils::Matrix2Hist(A, "Shaw_A",0.,1.,0.,1.);
  t.xTruth = MatrixUtils::Vec2Hist(x, 0., 1., "Shaw_x","Truth x fn.");
  t.bIdeal = MatrixUtils::Vec2Hist(b_ideal, 0., 1., "Shaw_b_ideal","Meas. b fn.");
  t.bNoisy = MatrixUtils::Vec2Hist(b, 0., 1., "Shaw_b","Meas. b fn.");
  for (int j=0; j<n; j++)
    t.bNoisy->SetBinError(j+1, noise);
  return t;
}

void
ShawSystem(const int n, TMatrixD &A, TVectorD &x, TVectorD &b,
           double noise)
{
  // Create a complete test problem for unfolding exercises that is
  // mathematically identical to shaw.m in the Matlab "regularization
  // tools" examples.
  //
  // To set up the test problem, just do, e.g.
  // TMatrixD A(n,n);
  // TVectorD xt(n), bt(n);      // true x and b (no noise)
  // ShawSystem(n, A, xt, bt);
  // Where n is an even integer.

  if (n%2)
  {
    gROOT->Error("ShawSystem()", "Even binning required");
    return;
  }

  // For computing A
  double h = TMath::Pi()/n;
  double co[n], psi[n];

  // For computing x vector
  double a1=2, a2=1, c1=6, c2=2, t1=0.8, t2=-0.5;

  for (int i=0; i<n; i++)
  {
    co[i] = TMath::Cos(-TMath::PiOver2() + (0.5 + i)*h);
    psi[i] = TMath::Pi()*TMath::Sin(-TMath::PiOver2() + (0.5 + i)*h);

    double arg1 = -c1*TMath::Power(-TMath::PiOver2() + (0.5 + i)*h - t1, 2);
    double arg2 = -c2*TMath::Power(-TMath::PiOver2() + (0.5 + i)*h - t2, 2);
    x(i) = a1*TMath::Exp(arg1) + a2*TMath::Exp(arg2);
  }

  for (int i=0; i<n/2; i++)
  {
    for (int j=i; j<n-i; j++)
    {
      double ss = psi[i] + psi[j];
      A(i,j) = TMath::Power((co[i]+co[j])*TMath::Sin(ss)/ss, 2);
      A(n-j-1,n-i-1) = A(i,j);
    }
    A(i,n-i-1) = TMath::Power(2*co[i], 2);
  }

  // Create upper triangular matrix from A
  TMatrixD Atri(A);
  for (int i=0; i<n; i++)
  {
    for (int j=i; j<n; j++)
    {
      if (j<=i)
        Atri(i,j)=0.;
    }
  }
  Atri.T();

  A += Atri;
  A *= h;

  b = A*x;

  // Add Gaussian white noise to b
  if (noise > 0.)
  {
    TRandom3 r3;
    for (int j=0; j<n; j++) b(j) += noise*r3.Gaus();
  }

  return;
}
