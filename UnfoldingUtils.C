
#include "UnfoldingUtils.h"
#include "MatrixUtils.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TF1.h"
#include <iostream>

using std::cout;
using std::endl;
using namespace MatrixUtils;

ClassImp(UnfoldingUtils);

UnfoldingUtils::UnfoldingUtils() :
  fVerbosity(1),
  fM(0),
  fN(0),
  fMeasX1(0.),
  fMeasX2(0.),
  fTrueX1(0.),
  fTrueX2(0.),
  fRegWeight(0.),
  fRegType(kNoReg),
  fHistAProb(0),
  fHistA(0),
  fHistATilde(0),
  fHistMeas(0),
  fHistMeasCov(0),
  fHistbTilde(0),
  fHistXini(0),
  fHistXtrue(0),
  fHistEff(0),
  fHistPrior(0)
{
}

UnfoldingUtils::UnfoldingUtils(TH2D *hA,
                               TH1D *hMeas,
                               TH2D *hMeasCov,
                               TH1D *hXini,
                               TH1D *hXtrue,
                               TH1D *hEff) :
  fVerbosity(1),
  fM(0),
  fN(0),
  fMeasX1(0.),
  fMeasX2(0.),
  fTrueX1(0.),
  fTrueX2(0.),
  fRegWeight(0.),
  fRegType(kNoReg),
  fHistAProb(0),
  fHistA(hA),
  fHistATilde(0),
  fHistMeas(hMeas),
  fHistMeasCov(hMeasCov),
  fHistbTilde(0),
  fHistXini(hXini),
  fHistXtrue(hXtrue),
  fHistEff(hEff),
  fHistPrior(0)
{
  fM = fHistA->GetNbinsX();    // Measured bins (rows)
  fN = fHistA->GetNbinsY();    // True/gen bins (cols)
  fTrueX1 = fHistA->GetYaxis()->GetXmin();
  fTrueX2 = fHistA->GetYaxis()->GetXmax();
  fMeasX1 = fHistA->GetXaxis()->GetXmin();
  fMeasX2 = fHistA->GetXaxis()->GetXmax();

  if (!BinningOk())
  {
    Error("UnfoldingUtils::UnfoldingUtils()", "Binning problem");
    gSystem->Exit(-1);
  }

  ComputeRescaledSystem();
}

void
UnfoldingUtils::ComputeRescaledSystem()
{
  // Matrix and vector members
  fMatA.ResizeTo(fM,fN);
  fMatAhat.ResizeTo(fM,fN);
  fMatATilde.ResizeTo(fM,fN);
  fMatB.ResizeTo(fM,fM);
  fMatBinv.ResizeTo(fM,fM);
  fVecb.ResizeTo(fM);
  fVecbErr.ResizeTo(fM);
  fVecbTilde.ResizeTo(fM);
  fVecXini.ResizeTo(fN);
  fVecXtrue.ResizeTo(fN);
  fSmoothingWeight.ResizeTo(fN);

  Info("UnfoldingUtils::ComputeRescaledSystem()","Initializing...");
  fMatA = Hist2Matrix(fHistA);
  fVecb = Hist2Vec(fHistMeas);

  if (fHistXini)
  {
    fVecXini = Hist2Vec(fHistXini);
  }
  else
  {
    fVecXini = Ones(fN);
    fHistXini = Vec2Hist(fVecXini, fTrueX1, fTrueX2,
                         "fHistXini", "Default x^{ini} (1.0)");
  }
  if (fHistXtrue)
    fVecXtrue = Hist2Vec(fHistXtrue);

  if (!fHistEff)
    fHistEff = Vec2Hist(Ones(fN), fTrueX1, fTrueX2,
                        "fHistEff", "Default efficiency (1.0)");
  if (!fHistPrior)
    fHistPrior = Vec2Hist(Ones(fN), fTrueX1, fTrueX2,
                          "fHistPrior", "Default prior (1.0)");
  // Probability matrix Ahat
  if (!fHistAProb)
  {
    fHistAProb = (TH2D *) fHistA->Clone("fHistAProb");
    NormalizeXSum(fHistAProb, fHistEff);
    fHistAProb->SetTitle("Probability matrix #hat{A}");
  }
  fMatAhat = Hist2Matrix(fHistAProb);

  // Data uncertainty
  for (int i=0; i<fM; i++)
    fVecbErr(i) = fHistMeas->GetBinError(i+1);

  // Create rescaled (~) quantities
  if (fHistMeasCov)
  {

    // If b has nontrivial covariance
    fMatB = Hist2Matrix(fHistMeasCov);

    // Eigendecomp: B = QRQ' where R_{ij} = r_i^2 \delta_{ij}
    // See eq. (33) (NIM A 372 (1996) 469-481)
    TDecompSVD svd(fMatB);
    TMatrixD Q = svd.GetU();   // Q=U=V if B is symm. & pos.-semidef
    TVectorD R = svd.GetSig(); // R(i,i)
    TVectorD r(fM);
    for (int i=0; i<fM; i++)
      r(i) = TMath::Sqrt(R(i));

    fMatATilde = DivColsByVector(Q*fMatA, r);
    fVecbTilde = ElemDiv(Q*fVecb, r);

    Info("UnfoldingUtils::ComputeRescaledSystem()",
         "Inverting covariance matrix...");
    fMatBinv = MoorePenroseInverse(fMatB);
  }
  else
  {
    // The usual case - b has indep. errors
    for (int i=0; i<fM; i++)
    {
      double var = fVecbErr(i)*fVecbErr(i);
      fMatB(i,i) = var;
      fMatBinv(i,i) = (var > 0) ? 1./var : 0;
    }
    fMatATilde = DivColsByVector(fMatA, fVecbErr);
    fVecbTilde = ElemDiv(fVecb, fVecbErr);
    fHistMeasCov = Matrix2Hist(fMatB, "fHistMeasCov",
                               fMeasX1,fMeasX2,fMeasX1,fMeasX2);
  }

  fHistATilde = Matrix2Hist(fMatATilde, "fHistATilde",
                            fMeasX1,fMeasX2,fTrueX1,fTrueX2);
  fHistATilde->SetTitle("covariance-scaled matrix #tilde{A}");

  fHistbTilde = Vec2Hist(fVecbTilde, fMeasX1, fMeasX2,
                         "fHistbTilde", "Scaled measured distribution");

  for (int i=0; i<fM; i++)
    fHistbTilde->SetBinError(i+1, 1.0);

  fSmoothingWeight = Ones(fN);

  Info("UnfoldingUtils::ComputeRescaledSystem()","Finished init step.");
  return;
}

bool
UnfoldingUtils::BinningOk()
{
  // Check for binning incompatibilities
  bool isok = true;

  int nMeas = fHistMeas->GetNbinsX(); // Must equal fM
  if (nMeas != fM)
  {
    gROOT->Warning("UnfoldingUtils::BinningOk()",
                   "Meas. bin mismatch: hMeas %d, TH2 (x-axis) %d",
                   nMeas, fM);
    isok = false;
  }
  if (fHistXini)
  {
    int nXini = fHistXini->GetNbinsX(); // Must equal fN
    if (nXini != fN)
    {
      gROOT->Warning("UnfoldingUtils::BinningOk()",
                     "True bin mismatch: x^ini %d, TH2 (y-axis) %d",
                     nXini, fN);
      isok = false;
    }
  }
  if (fHistXtrue)
  {
    int nXtrue = fHistXtrue->GetNbinsX(); // Must equal fN
    if (nXtrue != fN)
    {
      gROOT->Warning("UnfoldingUtils::BinningOk()",
                     "True bin mismatch: hXtrue %d, TH2 (y-axis) %d",
                     nXtrue, fN);
      isok = false;
    }
  }
  return isok;
}

TString
UnfoldingUtils::Algorithm(int type)
{
  TString s;
  switch (type)
  {
  case kSVDAlgo:       s = "SVD";      break;
  case kGSVDAlgo:      s = "GSVD";     break;
  case kRichLucyAlgo:  s = "RichLucy"; break;
  case kChi2MinAlgo:   s = "Chi2Min";  break;
  case kPCGLSAlgo:     s = "PCGLS";    break;
  default: s = "unrecognized_algorithm";
  }
  return s;
}

void
UnfoldingUtils::SetTrueRange(double x1, double x2)
{
  fTrueX1 = x1;
  fTrueX2 = x2;
}

void
UnfoldingUtils::SetPrior(TH1D *h)
{
  fHistPrior = (TH1D *)h->Clone("fHistPrior");
}

void
UnfoldingUtils::SetMeasRange(double x1, double x2)
{
  fMeasX1 = x1;
  fMeasX2 = x2;
}

TMatrixD
UnfoldingUtils::GetA(TString opt)
{
  if (opt.Contains("^"))
    return fMatAhat;
  else if (opt.Contains("~"))
    return fMatATilde;
  else
    return fMatA;
}

TVectorD
UnfoldingUtils::Getb(TString opt)
{
  if (opt.Contains("~"))
    return fVecbTilde;
  else
    return fVecb;
}

TH1D *
UnfoldingUtils::GetXTrueHist()
{
  if (!fHistXtrue)
    Error("UnfoldingUtils::GetXTrueHist()",
          "No xtrue histogram assigned for this problem");
  return fHistXtrue;
}

Double_t
UnfoldingUtils::GetbErrNorm()
{
  return TMath::Sqrt(fVecbErr*fVecbErr);
}

Double_t
UnfoldingUtils::GetbErrMean()
{
  return fVecbErr.Sum()/fVecbErr.GetNrows();
}

Double_t
UnfoldingUtils::GetbErrRMS()
{
  TVectorD e2 = ElemMult(fVecbErr,fVecbErr);
  return TMath::Sqrt(e2.Sum()/e2.GetNrows());
}


SVDResult
UnfoldingUtils::SVDAnalysis(TH2 *hA, TH1 *hb, TString opt)
{
  // Decompose A as U*Sigma*V' and study the components. If A is m x
  // n, then U is a column-orthogonal m x n matrix, Sigma is a
  // diagonal n x n matrix (stored as a vector), and V is
  // a column-orthonormal n x n matrix (V'*V = 1).

  static int id = 0; id++;
  SVDResult result;

  fTilde = (opt.Contains("~")) ? true : false;

  // Use stored members:
  // Select A, \hat{A}, or \tilde{A} using opt
  TMatrixD A = GetA(opt);
  TVectorD b = Getb(opt);

  // Or, use passed-in histograms
  if (hA)
  {
    A.ResizeTo(hA->GetNbinsX(), hA->GetNbinsY());
    A = Hist2Matrix(hA);
  }
  if (hb)
  {
    b.ResizeTo(hb->GetNbinsX());
    b = Hist2Vec(hb);
  }
  if (b.GetNrows()==0)
  {
    if (A.GetNrows())
      b.ResizeTo(A.GetNrows());
    else
      Warning("UnfoldingUtils::SVDAnalysis()",
              "Unspecified dimension of b vector");
  }

  TDecompSVD decomp(A);
  TVectorD sig = decomp.GetSig();
  TMatrixD U   = decomp.GetU();
  TMatrixD UT(TMatrixD::kTransposed, decomp.GetU());
  int ns = sig.GetNoElements();
  TVectorD utb = UT*b;

  if (utb.GetNrows() != sig.GetNrows())
    utb.ResizeTo(sig);

  for (int i=0; i<ns; i++)
    utb(i) = TMath::Abs(utb(i));
  TVectorD svc = ElemDiv(utb, sig); // SVD coefficients (abs. value)

  result.sigma = Vec2Hist(sig, 0., ns, Form("sig%d",id),
                          "#sigma_{i} ");
  result.UTb   = Vec2Hist(utb, 0., ns, Form("utb%d",id),
                          "#||{u^{T}_{i}#upointb} ");
  result.coeff = Vec2Hist(svc, 0., ns, Form("svc%d",id),
                          "#||{u^{T}_{i}#upointb} / #sigma_{i} ");
  result.U     = Matrix2Hist(U, Form("U_svd_%d",id),
                             0, ns, 0, ns);

  SetTH1Props(result.sigma, kBlack, 0, kBlack, kFullCircle, 1.0);
  SetTH1Props(result.UTb,   kBlue, 0, kBlue, kFullSquare, 1.0);
  SetTH1Props(result.coeff, kRed, 0, kRed, kOpenSquare, 1.0);

  return result;
}

GSVDResult *
UnfoldingUtils::GSVDAnalysis(TMatrixD &L, double lambda, TH2 *hA, TH1 *hb,
                             TString opt)
{
  // Decompose A, L jointly as A = U*C*X', L = V*S*X'.
  // If A is m x n, and L is p x n, and A has full rank,
  // U  is m x n
  // V  is p x p
  // X' is n x n
  // C  is n x n
  // S  is p x n
  // alpha and beta are the "interesting" diagonal elements of C,S
  // respectively, and have length p.

  static int id = 0; id++;
  GSVDResult *gsvd = new GSVDResult();

  // Use stored members:
  // Select A, \hat{A}, or \tilde{A} ("", "^", or "~")
  TMatrixD A = GetA(opt);
  TVectorD b = Getb(opt);

  // Or, use passed-in histograms (supercedes opt)
  if (hA)
  {
    A.ResizeTo(hA->GetNbinsX(), hA->GetNbinsY());
    A = Hist2Matrix(hA);
  }
  if (hb)
  {
    b.ResizeTo(hb->GetNbinsX());
    b = Hist2Vec(hb);
  }

  int m = A.GetNrows();
  int n = A.GetNcols();
  int p = L.GetNrows();

  GSVDecompResult g = GSVD(A,L);
  TMatrixD UT(TMatrixD::kTransposed, g.U);
  TMatrixD X(TMatrixD::kInverted, g.XT);
  TVectorD utb(UT*b);

  // GSVD expansion coeffs
  TVectorD c = ElemDiv(utb,g.alpha);

  // Tikhonov filter factors
  TVectorD f(n);
  for (int i=0; i<n; i++)
  {

    if (i >= n-p)
    {
      double g2 = g.gamma(i)*g.gamma(i);
      f(i) = g2 / (g2 + lambda*lambda);
    }
    else
      f(i) = 1.0;
  }

  // Compute F*C^+ as a step toward computing the regularized inverse
  // of A. Since C is diagonal, save time by computing its
  // pseudoinverse without using SVD.
  TMatrixD FCd(n,n);
  for (int i=0; i<n; i++)
  {
    double val = g.C(i,i);
    FCd(i,i) = (val > 0.) ? f(i)/val : 0.;
  }
  gsvd->Ap.ResizeTo(n,m);
  gsvd->Ap = X*FCd*UT;

  // Covariance matrices:
  // Cov(w) = Ap*B*Ap'
  // Cov(x)_ik = xini_i Cov(w)_ik xini_k
  gsvd->covw.ResizeTo(n,n);
  gsvd->covx.ResizeTo(n,n);
  gsvd->covxInv.ResizeTo(n,n);

  TMatrixD B(fMatB);
  if (opt.Contains("~"))
    B.UnitMatrix();
  gsvd->covb.ResizeTo(m,m);
  gsvd->covb = B;

  // Compute cov(x)
  TMatrixD tmp(B, TMatrixD::kMultTranspose, gsvd->Ap);
  gsvd->covw = gsvd->Ap * tmp;
  for (int i=0; i<n; i++)
  {
    for (int k=0; k<n; k++)
    {
      gsvd->covx(i,k) = fVecXini(i) * gsvd->covw(i,k) * fVecXini(k);
    }
  }

  // Compute cov(x)^-1 (without using cov(x), see Hocker eq. 53)
  double AAjk = 0;
  for (int j=0; j<n; j++)
  {
    for (int k=0; k<n; k++)
    {
      AAjk = 0;
      for (int i=0; i<n; i++)
      {
        AAjk += A(i,j)*A(i,k);
      }
      if (fVecXini(j)==0 || fVecXini(k)==0)
        gsvd->covxInv(j,k) = 0;
      else
        gsvd->covxInv(j,k) = AAjk/fVecXini(j)/fVecXini(k);
    }
  }

  TVectorD regc = ElemMult(f,c);
  TVectorD wreg = X*regc;
  TVectorD xreg = ElemMult(fVecXini, wreg);

  // Save to output for parameter optimization analysis
  gsvd->bInc.ResizeTo(m);
  gsvd->bInc = (LMatrix(m,kUnitMatrix) - g.U * UT)*b;

  // Assign output struct members
  gsvd->n      = n;
  gsvd->m      = m;
  gsvd->p      = p;
  gsvd->lambda = lambda;
  gsvd->alpha.ResizeTo(n);  gsvd->alpha  = g.alpha;
  gsvd->beta .ResizeTo(n);  gsvd->beta   = g.beta;
  gsvd->gamma.ResizeTo(n);  gsvd->gamma  = g.gamma;
  gsvd->f    .ResizeTo(n);  gsvd->f      = f;
  gsvd->UTb  .ResizeTo(n);  gsvd->UTb    = utb;
  gsvd->coeff.ResizeTo(n);  gsvd->coeff  = c;
  gsvd->regc .ResizeTo(n);  gsvd->regc   = regc;

  // Copy A,L,and b to output to ensure that the correct quantities
  // are used later.
  gsvd->X.ResizeTo(X);      gsvd->X = X;
  gsvd->U.ResizeTo(g.U);    gsvd->U = g.U;
  gsvd->V.ResizeTo(g.V);    gsvd->V = g.V;
  gsvd->L.ResizeTo(L);      gsvd->L = L;
  gsvd->A.ResizeTo(A);      gsvd->A = A;
  gsvd->b.ResizeTo(b);      gsvd->b = b;

  gsvd->UHist  = Matrix2Hist(g.U, Form("U_gsvd_%d",id),
                             fTrueX1, fTrueX2,0,n);
  gsvd->XHist  = Matrix2Hist(X, Form("X_gsvd_%d",id),
                             fTrueX1, fTrueX2,0,n);
  gsvd->wregHist = Vec2Hist(wreg, fTrueX1,fTrueX2,
                            Form("gsvd_wreg_%d",id),
                            Form("w (#lambda = %g)", lambda));
  gsvd->xregHist = Vec2Hist(xreg, fTrueX1,fTrueX2,
                            Form("gsvd_xreg_%d",id),
                            Form("x (#lambda = %g)", lambda));

  gsvd->bregHist = Vec2Hist(fMatAhat*xreg, fMeasX1,fMeasX2,
                            Form("gsvd_breg_%d",id),
                            Form("Ax_{#lambda} (#lambda = %g)", lambda));

  // Assign uncertainties
  for (int j=0; j<n; j++)
  {
    gsvd->wregHist->SetBinError(j+1, TMath::Sqrt(gsvd->covw(j,j)));
    gsvd->xregHist->SetBinError(j+1, TMath::Sqrt(gsvd->covx(j,j)));
  }

  // Absolute values for plotting
  TVectorD utbAbs(utb);
  TVectorD cAbs(c);
  TVectorD rcAbs(ElemMult(f,utb));   // TVectorD rcAbs(regc);
  for (int i=0; i<n; i++)
  {
    if (utbAbs(i) < 0) utbAbs(i) *= -1;
    if (cAbs(i) < 0)     cAbs(i) *= -1;
    if (rcAbs(i) < 0)   rcAbs(i) *= -1;
  }
  gsvd->UTbAbs   = Vec2Hist(utbAbs, 0,n,Form("gsvd_utb_abs%d",id),
                            "#||{u^{T}#upointb} ");
  gsvd->coeffAbs = Vec2Hist(cAbs, 0,n,Form("gsvd_c_abs%d",id),
                            "#||{u^{T}#upointb}/#alpha ");
  gsvd->regcAbs  = Vec2Hist(rcAbs, 0,n,Form("gsvd_rc_abs%d",id),
                            "f#||{u^{T}#upointb} ");

  SetTH1Props(gsvd->UTbAbs,   kBlue, 0, kBlue, kFullSquare, 1.0);
  SetTH1Props(gsvd->coeffAbs, kRed, 0, kRed, kOpenSquare, 1.0);
  SetTH1Props(gsvd->regcAbs, kMagenta+1, 0, kMagenta+1, kOpenCircle, 1.0);
  SetTH1Props(gsvd->xregHist, kGreen+2, 0, kGreen+2, kFullCircle, 1.0);

  return gsvd;
}

TCanvas *
UnfoldingUtils::DrawSVDPlot(SVDResult svdhists, double ymin, double ymax,
                            TString opt)
{
  static int i=0; i++;
  TCanvas *c = new TCanvas(Form("csvd%d",i), Form("csvd%d",i), 1);

  // Draw frame histogram to set limits, title, etc.
  int nx = svdhists.sigma->GetNbinsX();
  TH1F *hsvd = new TH1F(Form("hsvd%d",i), "SV Components;column index i;",
                        200, 0, nx);
  hsvd->Draw();
  hsvd->GetYaxis()->SetRangeUser(ymin, ymax);

  // Singular value spectrum
  if (opt.Contains("sig"))
    svdhists.sigma->Draw("plsame");

  // Draw |U'*b| and |U'*b|/sigma
  svdhists.UTb->Draw("plsame");
  svdhists.coeff->Draw("plsame");
  gPad->SetLogy();

  TLegend *leg = new TLegend(0.75, 0.75, 0.99, 0.99);
  if (opt.Contains("sig"))
    leg->AddEntry(svdhists.sigma, svdhists.sigma->GetTitle(), "p");
  leg->AddEntry(svdhists.UTb, svdhists.UTb->GetTitle(), "ep");
  leg->AddEntry(svdhists.coeff, svdhists.coeff->GetTitle(), "ep");
  leg->SetFillColor(kNone);
  leg->Draw();

  return c;
}

TCanvas *
UnfoldingUtils::DrawGSVDPlot(GSVDResult *gsvd, double ymin, double ymax,
                             TString opt)
{
  static int i=0; i++;
  TCanvas *c = new TCanvas(Form("cgsvd%d",i), Form("cgsvd%d",i), 1);
  TLegend *leg = new TLegend(0.75, 0.75, 0.99, 0.99,
                             Form("NDF_{eff} = %.1f",
                                  gsvd->f.Sum()));
  int nx = gsvd->UTbAbs->GetNbinsX();
  TH1F *h = new TH1F(Form("hgsvd%d",i),
                     Form("GSVD Components (#lambda = %g);column index i;",
                          gsvd->lambda), 200, 0, nx);
  h->Draw();
  h->GetYaxis()->SetRangeUser(ymin, ymax);
  gPad->SetLogy();
  gsvd->UTbAbs->Draw("plsame");
  gsvd->regcAbs->Draw("plsame");


  leg->AddEntry(gsvd->UTbAbs, gsvd->UTbAbs->GetTitle(), "ep");
  leg->AddEntry(gsvd->regcAbs, gsvd->regcAbs->GetTitle(), "ep");

  // Draw coeffs also if requested
  if (opt.Contains("c"))
  {
    gsvd->coeffAbs->Draw("plsame");
    leg->AddEntry(gsvd->coeffAbs, gsvd->coeffAbs->GetTitle(), "ep");
    // TODO: add damped coeffs f|U'*b|/alpha
  }

  leg->SetFillColor(kNone);
  leg->Draw();

  if (opt.Contains("~"))
  {
    TLine one;
    TAxis *ax = h->GetXaxis();
    one.DrawLine(ax->GetXmin(), 1.0, ax->GetXmax(), 1.0);
  }

  return c;
}

double
UnfoldingUtils::SmoothingNorm(TVectorD &x, int regtype)
{
  double sn = 0;
  switch (regtype)
  {
  case kNone:
    sn = 0;
    break;
  case kTotCurv:
    sn = Curvature(x);
    break;
  }
  return sn;
}

double
UnfoldingUtils::Curvature(TVectorD &x)
{
  // Eq. (38), NIM A 372 (1996) 469-481.
  double delta=0, val=0;
  for (int i=1; i<fN-1; i++)
  {
    delta = fSmoothingWeight(i)*(x(i+1) - x(i)) - (x(i) - x(i-1));
    val += delta*delta;
  }
  return val;
}

double
UnfoldingUtils::RegChi2(const double *pars)
{
  // Returns a modified chi squared value. Designed to be called by
  // TMinuit for minimization of the return value.

  double chi2 = 0;
  // Fit parameters vector
  TVectorD x(fN);
  for (int i=0; i<fN; i++)
    x(i) = pars[i];

  // Unmodified chi^2 (Ax-b)'*Binv*(Ax-b)
  if (fTilde)
  {
    TVectorD resid = fMatATilde*x - fVecbTilde;
    chi2 = resid*resid;
  }
  else
  {
    TVectorD resid = fMatAhat*x - fVecb;
    chi2 = resid*(fMatBinv*resid);
  }

  // Additive chi^2 modifier (reg. penalty value)
  double mod = fRegWeight*SmoothingNorm(x, fRegType);

  chi2 += mod;
  return chi2;
}

UnfoldingResult
UnfoldingUtils::UnfoldChiSqMin(TVectorD &regWts, TString opt)
{
  // Unfold by minimizing return values of RegChi2().
  // To seed the fitter, call SetPrior(TH1*).

  static int id=0; id++;
  int nRegWts = regWts.GetNrows();
  int nPars = fN;
  UnfoldingResult result;

  Info("UnfoldingUtils::UnfoldChiSqMin()",
       "Using reg type %d, weight %g. Option: %s",
       fRegType, fRegWeight, opt.Data());

  // Control which A and b are used in RegChi2().
  if (opt.Contains("~"))
    fTilde = true;

  double kbins[nRegWts+1], xbins[fN+1];
  for (int j=0; j<=fN; j++)
  {
    xbins[j] = fTrueX1 + j*(fTrueX2-fTrueX1)/fN;
  }
  for (int k=0; k<nRegWts; k++)
  {
    kbins[k] = regWts(k);
  }
  kbins[nRegWts] = 2*kbins[nRegWts-1]-kbins[nRegWts-2];

  //  result.XRegHist = Matrix2Hist(result.XReg, XRegErr, Form("xregHist_%d",id), xbins, kbins);
  result.XRegHist = new TH2D(Form("hChsq%d",id), Form("Minimum #chi^{2}_{reg} solutions"),
                             fN, xbins, nRegWts, kbins);

  result.XRegHist->GetXaxis()->CenterTitle();
  result.XRegHist->GetYaxis()->CenterTitle();
  result.XRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.XRegHist->GetYaxis()->SetTitleOffset(1.8);
  result.XRegHist->SetTitle("GSVD solutions: x_{#lambda};x;#lambda");

  result.LCurve = new TGraph();
  result.LCurve->SetTitle("TMinuit L-Curve;Unregulated total #chi^{2};total curvature");

  // Set up the chi^2 minimizer
  //  TFitterMinuit* min = new FitterMinuit();
  ROOT::Math::Minimizer *min =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(this, &UnfoldingUtils::RegChi2, nPars);
  min->SetFunction(f);

  // Initialize tmx array and pass to minimizer.
  double *tmx = new double[fN];
  double stepSize = 0.1;

  for (int k=0; k<nRegWts; k++)
  {
    for (int j=0; j<fN; j++)
    {
      tmx[j] = fHistPrior->GetBinContent(j+1);

      // Require all parameters to have a minimum positive value
      if (tmx[j] < 0)
      {
        tmx[j] = 1e-6;
      }

      // Constrain pars to be >= 0 unless requested
      if (opt.Contains("-"))
        min->SetVariable(j, Form("xpar%d",j), tmx[j], stepSize);
      else
        min->SetLowerLimitedVariable(j, Form("xpar%d",j), tmx[j], stepSize, 0.0);

      if (0) Printf("%g", tmx[j]);
    }

    fRegWeight = regWts(k);

    Info("UnfoldingUtils::UnfoldChiSqMin()",
         "Initial (regularized) chi squared = %g for reg. weight %g",
         RegChi2(tmx), regWts(k));

    min->Minimize();

    TVectorD w(fN); // Solution for this k
    for (int j=0; j<nPars; j++)
    {
      tmx[j] = min->X()[j];
      double val = tmx[j];
      double err = val < 1 ? 1. : min->Errors()[j];
      w(j) = val;

      result.XRegHist->SetBinContent(j+1,k+1,val);
      result.XRegHist->SetBinError(j+1,k+1,err);

    }
    fRegWeight = 0.;
    result.LCurve->SetPoint(k, RegChi2(tmx), Curvature(w));
  } // k loop

  return result;
}

TMatrixD
UnfoldingUtils::RegularizedInverseResponse(GSVDResult *gsvd, double lambda)
{
  // Compute A^# from GSVD components as X * F_lambda * C^dagger * U'
  // where F is the diagonal matrix of filter factors for this lambda value.
  int n  = gsvd->n;
  int m  = gsvd->m;
  int p  = gsvd->p;
  TMatrixD FCd(n,n);      // Filter factor matrix * pseudoinverse(diag(alpha))
  TMatrixD Ap(n,m);       // Regularized inverse A^#
  TMatrixD UT(TMatrixD::kTransposed, gsvd->U);

  // Create Tikhonov filter factors for this lambda
  TVectorD f(n);
  for (int i=0; i<n; i++)
  {
    if (i >= n-p)
    {
      double g2 = gsvd->gamma(i)*gsvd->gamma(i);
      f(i) = g2 / (g2 + lambda*lambda);
    }
    else
      f(i) = 1.0;
  }
  for (int i=0; i<fN; i++)
  {
    FCd(i,i) = (gsvd->alpha(i) > 0.) ? f(i)/gsvd->alpha(i) : 0.;
  }
  Ap = gsvd->X * FCd * UT;
  return Ap;
}

UnfoldingResult
UnfoldingUtils::UnfoldTikhonovGSVD(GSVDResult *gsvd,
                                   TVectorD &lambda,
                                   TString /*opt*/)
{
  UnfoldingResult result;
  static int id=0; id++; // So this fn. can be called more than once

  int nk = lambda.GetNrows();
  int n  = gsvd->n;
  int m  = gsvd->m;
  int p  = gsvd->p;
  TMatrixD WRegErr(n,nk);
  TMatrixD XRegErr(n,nk);
  result.WReg.ResizeTo(n, nk);
  result.XReg.ResizeTo(n, nk);

  result.LCurve    = new TGraph(nk);
  result.GcvCurve  = new TGraph(nk);
  result.RhoCurve  = new TGraph(nk);
  result.FilterSum = new TGraph(nk);
  TGraph *fsInv    = new TGraph(nk); // FilterSum w/ x,y swapped


  result.hwCov = new TH3D(Form("hwCov_gsvd_%d", id),
                          Form("hwCov_gsvd_%d", id),
                          fN,fTrueX1,fTrueX2,
                          fN,fTrueX1,fTrueX2,
                          nk, lambda(0), lambda(nk-1));
  result.hxCov = new TH3D(Form("hxCov_gsvd_%d", id),
                          Form("hxCov_gsvd_%d", id),
                          fN,fTrueX1,fTrueX2,
                          fN,fTrueX1,fTrueX2,
                          nk, lambda(0), lambda(nk-1));

  result.lambdaGcv = 0;
  double gcvMin = 1e99;
  result.lambdaRho = 0;
  double rhoMin = 1e99;
  result.lambdaLcv = 0;
  result.kStf = 0;
  result.lambdaStf = 0;

  // Stuff for computing covariance
  TMatrixD B(gsvd->covb); // Error matrix of b
  TMatrixD FCd(n,n);      // Filter factor matrix * pseudoinverse(diag(alpha))
  // TMatrixD Ap(n,m);       // Regularized inverse A^#
  TMatrixD UT(TMatrixD::kTransposed, gsvd->U);

  // Compute the number of significant GSVD coefficients
  int nSignificantGSVDCoeffs = 0;
  double errThreshold = 0;
  for (int i=0; i<fM; i++)
    errThreshold += B(i,i)/fM;

  for (int j=0; j<fN-1; j++)
  {
    if (gsvd->UTbAbs->GetBinContent(j+1) > errThreshold &&
        gsvd->UTbAbs->GetBinContent(j+2) < errThreshold)
    {
      nSignificantGSVDCoeffs = j;
      break;
    }
  }
  Printf("# coeffs > %.2f = %d", errThreshold, nSignificantGSVDCoeffs);

  // Scan over lambda values, generate nk solutions
  for (int k=0; k<nk; k++)
  {

    // Create Tikhonov filter factors for this lambda
    double l = lambda(k);
    TVectorD f(n);
    for (int i=0; i<n; i++)
    {
      if (i >= n-p)
      {
        double g2 = gsvd->gamma(i)*gsvd->gamma(i);
        f(i) = g2 / (g2 + l*l);
      }
      else
        f(i) = 1.0;
    }

    // Damped GSVD coefficients
    TVectorD regc = ElemMult(f, gsvd->coeff);
    TVectorD wreg = gsvd->X*regc;
    TVectorD xreg = ElemMult(fVecXini, wreg);

    // Regularized Inverse A^# for this lambda
    TMatrixD Ap = RegularizedInverseResponse(gsvd, lambda(k));

    // Covariance matrices for this lambda:
    // Cov(w) = Ap*B*Ap'
    // Cov(x)_ik = xini_i Cov(w)_ik xini_k
    TMatrixD covw = Ap * TMatrixD(B, TMatrixD::kMultTranspose, Ap);
    TMatrixD covx(covw);

    // For each lambda, cov(w) and cov(x) are added to TH3s.
    for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
      {
        double cw = covw(i,j);
        double cx = fVecXini(i) * cw * fVecXini(j);
        covx(i,j) = cx;
        result.hwCov->SetBinContent(i+1,j+1,k+1,cw);
        result.hxCov->SetBinContent(i+1,j+1,k+1,cx);
        if (i==j) WRegErr(j,k) = TMath::Sqrt(cw);
        if (i==j) XRegErr(j,k) = TMath::Sqrt(cx);
      }
    }

    // Regularized solution
    for (int j=0; j<n; j++)
    {
      result.WReg(j,k) = wreg(j);
      result.XReg(j,k) = xreg(j);
    }

    // Parameter optimization analysis -------------------------------
    // ---------------------------------------------------------------

    // Compute Lx_reg
    TVectorD vec = ElemDiv(gsvd->UTb, gsvd->gamma);
    vec = ElemMult(f, vec);
    TVectorD Lx = vec.GetSub(n-p, n-1);
    Lx = gsvd->V * Lx;
    double lxnorm = TMath::Sqrt(Lx*Lx);

    // Compute r = b - Ax_reg
    TVectorD f1(f);
    for (int j=0; j<n; j++)
      f1(j) = 1-f(j);
    vec = ElemMult(f1, gsvd->UTb);

    TMatrixD Up = gsvd->U.GetSub(0,m-1,n-p,n-1); // m x n --> m x p
    TVectorD r = Up * vec.GetSub(n-p, n-1) - gsvd->bInc;
    double rnorm = TMath::Sqrt(r*r);
    double rho = 0;
    double fsum = f.Sum();
    double gcv = rnorm / (m - fsum);

    // Compute mean global correlation coefficients (V. Blobel)
    //    double rho_j = 0;
    double rhoMean = 0;
    double rhoj2 = 0;
    TMatrixD covxInv = MoorePenroseInverse(covx);

    for (int j=0; j<n; j++)
    {
      double vxreg = covx(j,j);
      double vxinv = covxInv(j,j);
      //      double vxinv = (gsvd->covxInv)(j,j);
      double prod = vxreg*vxinv;

      // if (prod == 0 || prod < 1.)
      //  rho_j = 0;
      // else
      //  rho_j = TMath::Sqrt(1 - 1./prod);

      // rhoMean += rho_j / n;

      if (prod == 0)
        rhoj2 = 0;
      else
        rhoj2 = 1 - 1./prod;

      rhoMean += rhoj2 / n;
    }
    rhoMean = TMath::Sqrt(rhoMean);
    rho = rhoMean;

    result.LCurve->SetPoint(k, rnorm, lxnorm);
    result.GcvCurve->SetPoint(k, lambda(k), gcv);
    result.RhoCurve->SetPoint(k, lambda(k), rhoMean);
    result.FilterSum->SetPoint(k, lambda(k), fsum);
    fsInv->SetPoint(k, fsum, lambda(k));

    // if (fsum > nSignificantGSVDCoeffs) {
    //   result.lambdaStf = lambda(k);
    //   result.kStf = k;
    // }
    if (gcv < gcvMin)
    {
      gcvMin = gcv;
      result.lambdaGcv = lambda(k);
      result.kGcv = k;
    }
    if (rho < rhoMin)
    {
      rhoMin = rho;
      result.lambdaRho = lambda(k);
      result.kRho = k;
    }
  }

  result.lambdaStf = fsInv->Eval(double(nSignificantGSVDCoeffs));
  result.kStf = TMath::BinarySearch(nk, result.FilterSum->GetX(),
                                    result.lambdaStf);

  // Assign result.Lcurvature and result.kLcv
  result.LCurvature = LogCurvature(result.LCurve, lambda, result.kLcv);
  result.lambdaLcv = lambda(result.kLcv);

  int kg = result.kGcv+1;
  result.RhoCurve->SetName("gsvd_rho");
  result.RhoCurve->SetTitle("GSVD mean global correlation coefficients;"
                            "#lambda;#LT#rho#GT");
  result.LCurve->SetTitle("GSVD L-Curve;||Ax_{#lambda}-b||_{2};"
                          "||Lx_{#lambda}||_{2}");

  result.LCurvature->SetName("gsvd_lcc");
  result.LCurvature->SetTitle("GSVD L-Curve log curvature;#lambda;"
                              "log curvature");
  result.GcvCurve->SetNameTitle("gsvd_gcv","GSVD cross-validation curve;"
                                "#lambda;G(#lambda)");
  result.FilterSum->SetNameTitle("gsvd_fsum", "GSVD Tikhonov filter factor sum"
                                 " (= effective NDF);#lambda;effective NDF");

  result.RhoCurve->GetXaxis()->CenterTitle();
  result.RhoCurve->GetYaxis()->CenterTitle();
  result.RhoCurve->GetXaxis()->SetTitleOffset(1.3);
  result.RhoCurve->GetYaxis()->SetTitleOffset(1.3);

  result.LCurve->GetXaxis()->CenterTitle();
  result.LCurve->GetYaxis()->CenterTitle();
  result.LCurve->GetXaxis()->SetTitleOffset(1.3);
  result.LCurve->GetYaxis()->SetTitleOffset(1.3);

  result.LCurvature->GetXaxis()->CenterTitle();
  result.LCurvature->GetYaxis()->CenterTitle();
  result.LCurvature->GetXaxis()->SetTitleOffset(1.3);
  result.LCurvature->GetYaxis()->SetTitleOffset(1.3);

  result.GcvCurve->GetXaxis()->CenterTitle();
  result.GcvCurve->GetYaxis()->CenterTitle();
  result.GcvCurve->GetXaxis()->SetTitleOffset(1.3);
  result.GcvCurve->GetYaxis()->SetTitleOffset(1.3);

  result.FilterSum->GetXaxis()->CenterTitle();
  result.FilterSum->GetYaxis()->CenterTitle();
  result.FilterSum->GetXaxis()->SetTitleOffset(1.3);
  result.FilterSum->GetYaxis()->SetTitleOffset(1.3);

  double kbins[nk+1], xbins[fN+1];
  for (int j=0; j<=fN; j++)
  {
    xbins[j] = fTrueX1 + j*(fTrueX2-fTrueX1)/fN;
  }
  for (int k=0; k<nk; k++)
  {
    kbins[k] = lambda(k);
  }
  kbins[nk] = 2*kbins[nk-1]-kbins[nk-2];

  result.WRegHist = Matrix2Hist(result.WReg, WRegErr, 
                                Form("wregHist_%d",id), xbins, kbins);
  result.WRegHist->GetXaxis()->CenterTitle();
  result.WRegHist->GetYaxis()->CenterTitle();
  result.WRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.WRegHist->GetYaxis()->SetTitleOffset(1.8);
  result.WRegHist->SetTitle("GSVD solutions: w_{#lambda};w;#lambda");


  result.XRegHist = Matrix2Hist(result.XReg, XRegErr, 
                                Form("xregHist_%d",id), xbins, kbins);
  result.XRegHist->GetXaxis()->CenterTitle();
  result.XRegHist->GetYaxis()->CenterTitle();
  result.XRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.XRegHist->GetYaxis()->SetTitleOffset(1.8);
  result.XRegHist->SetTitle("GSVD solutions: x_{#lambda};x;#lambda");

  int kr = result.kRho+1;
  int kl = result.kLcv+1;
  int ks = result.kStf+1;

  result.hGcv =
    result.XRegHist->ProjectionX(Form("gsvd_%d_gcv_bin%d",id,kg),kg,kg);

  result.hRho =
    result.XRegHist->ProjectionX(Form("gsvd_%d_rho_bin%d",id,kr),kr,kr);

  result.hLcv =
    result.XRegHist->ProjectionX(Form("gsvd_%d_lcv_bin%d",id,kl),kl,kl);

  result.hStf =
    result.XRegHist->ProjectionX(Form("gsvd_%d_stf_bin%d",id,ks),ks,ks);

  return result;
}

TH1D *
UnfoldingUtils::UnfoldSVD(double lambda,
                          TObjArray *output,
                          TString opt,
                          TH2 *hA,
                          TH1 *hb,
                          TH1 *hXini)
{
  static int id=0; id++; // So this fn. can be called more than once

  int matrixType = k2DerivBCR; // favor reflected w at boundaries
  if (opt.Contains("BC0"))
    matrixType = k2DerivBC0;   // favor w=0 at boundaries
  if (opt.Contains("I"))
    matrixType = kUnitMatrix;   // favor w=0 at boundaries

  fTilde = (opt.Contains("~")) ? true : false;

  // Setup
  TMatrixD A(fMatA);
  TVectorD b(fVecb);
  TVectorD xini(fVecXini); // All 1's if no fHistXini

  if (fTilde)
  {
    A = fMatATilde;
    b = fVecbTilde;
  }
  else
  {
    A = fMatAhat;
  }

  // Optionally, use passed-in histos instead of members. Useful for
  // toy MC trials in computing covariance matrix.
  if (hA)
    A = Hist2Matrix(hA);
  if (hb)
    b = Hist2Vec(hb);
  if (hXini)
    xini = Hist2Vec(hXini);

  TMatrixD L    = LMatrix(A.GetNcols(), matrixType, 1e-5);
  TMatrixD Linv = MoorePenroseInverse(L);
  TMatrixD LTi(L);
  LTi.T(); // Inverse transpose L^{-T}
  LTi = MoorePenroseInverse(LTi);

  // Compute SVD of AL^{-1} & get results
  TDecompSVD decomp(A*Linv);
  TVectorD s_vec   = decomp.GetSig();
  TMatrixD UT      = decomp.GetU(); UT.T();
  TMatrixD V       = decomp.GetV();
  TMatrixD VT(V); VT.T();
  TVectorD d       = UT*b;             // eq. 44
  int nd           = d.GetNrows();
  TVectorD tf(nd);   // Tikhonov filter factors
  TVectorD z(nd);    // z = d/s * tf
  TVectorD dz(nd);

  // Abs. value vectors - for analysis
  TVectorD absd(nd); // absd(i) is |(U^T*b)_i|
  TVectorD dlam(nd); // dlam(i) is tf_i*|d_i|, or s_i*|z_i|

  // Compute filter factors and z
  double s0 = s_vec(0);
  double si = 0;
  double smin = 1e-12*s0;
  TMatrixD Z(nd,nd);
  for (int i=0; i<nd; i++)
  {
    si = s_vec(i);
    if (si < smin) si = smin;
    tf(i) = si*si/(si*si + lambda*lambda);
    z(i) = d(i)/si * tf(i);
    dz(i) =  tf(i)/si;
    Z(i,i) = dz(i)*dz(i);

    // Extra - not part of solution
    absd(i) = TMath::Abs(d(i));
    dlam(i) = absd(i)*tf(i);
  }

  // Compute final solutions
  TVectorD w = Linv * V * z;
  TVectorD x(w);
  x = ElemMult(xini,w);

  TVectorD resid = A*x - b;

  // and covariance matrices:
  TMatrixD Wcov = Linv*V*Z*VT*LTi;
  TMatrixD Xcov(nd,nd);
  for (int i=0; i<nd; i++)
  {
    for (int k=0; k<nd; k++)
    {
      Xcov(i,k) = xini(i)*Wcov(i,k)*xini(k);
    }
  }

  // X inverse (without using Xcov, see eq. 53)
  TMatrixD Xinv(nd, nd);
  double AAjk = 0;
  for (int j=0; j<nd; j++)
  {
    for (int k=0; k<nd; k++)
    {
      AAjk = 0;
      for (int i=0; i<nd; i++)
      {
        AAjk += A(i,j)*A(i,k);
      }
      Xinv(j,k) = AAjk/xini(j)/xini(k);
    }
  }

  // Add components to output list
  if (output)
  {
    TH1D *hs = Vec2Hist(s_vec, 0., nd, Form("hs%d",id), "s_{i} ");
    TH1D *hd = Vec2Hist(absd,  0., nd, Form("hd%d",id), 
                        Form("#||{d_{i}} = #||{(U^{T}b)_{i}}"));
    TH1D *hl = Vec2Hist(dlam,  0., nd, Form("hl%d",id), 
                        Form("#||{d^{(#lambda)}_{i}}, #lambda = %g ",lambda));
    TH1D *hw = Vec2Hist(w,     0., nd, Form("hw%d",id), 
                        Form("w^{#lambda = %g} ",lambda));
    TH1D *ht = Vec2Hist(tf,    0., nd, Form("ht%d",id), 
                        Form("s_{i}^{2}/(s_{i}^{2}+#lambda^{2}), #lambda = %g ",
                             lambda));
    TH1D *hr = Vec2Hist(resid, 0., nd, Form("hr%d",id), 
                        Form("Ax^{(#lambda)}-b, #lambda = %g ",lambda));

    SetTH1Props(hs, kBlack, 0, kBlack, kFullCircle, 1.0);
    SetTH1Props(hd, kBlue,  0, kBlue,  kFullSquare, 1.0);
    SetTH1Props(hl, kMagenta+1,  0, kMagenta+1,  kOpenSquare, 1.0);
    SetTH1Props(hw, kCyan+2,0, kCyan+2,kOpenCircle, 1.0);
    SetTH1Props(ht, kBlack, 0, kBlack, kOpenCircle, 1.0);

    // Covariance matrix of solution and its inverse
    TH2D *hWcov = Matrix2Hist(Wcov, Form("hWcov%d",id), 0, nd, 0, nd);
    TH2D *hXcov = Matrix2Hist(Xcov, Form("hXcov%d",id), 0, nd, 0, nd);
    TH2D *hXinv = Matrix2Hist(Xinv, Form("hXinv%d",id), 0, nd, 0, nd);

    for (int i=0; i<nd; i++)
    {
      hd->SetBinError(i+1, 1.0);
      hl->SetBinError(i+1, 1.0);
    }
    output->Add(hs);
    output->Add(hd);
    output->Add(hl);
    output->Add(hw);
    output->Add(ht);
    output->Add(hr);
    output->Add(hWcov);
    output->Add(hXcov);
    output->Add(hXinv);
  }

  TH1D *hx = Vec2Hist(x,fTrueX1,fTrueX2,Form("hSVD%d",id),
                      Form("#lambda = %g", lambda));

  for (int i=0; i<nd; i++)
  {
    hx->SetBinError(i+1, TMath::Sqrt(Xcov(i,i)));
  }

  return hx;
}

TH1D *
UnfoldingUtils::UnfoldTLS()
{
  TMatrixD A = GetA();
  TVectorD b = Getb();

  TMatrixD C(fM, fN+1); // C = [A b]
  C.SetSub(0,0,A);
  for (int i=0; i<fM; i++) C(i,fN) = b(i);

  TDecompSVD svdc(C);
  TMatrixD V = svdc.GetV();
  TVectorD x(fN);
  for (int i=0; i<fN; i++) x(i) = V(i, fN);

  x *= -1./V(fN,fN);
  x = ElemMult(x, fVecXini);

  //  TMatrixD VAB = V.GetSub(0,fN-1,fN,fN);
  // TMatrixD VBB = V.GetSub(fN,fN,fN,fN);
  // TMatrixD VBBinv(TMatrixD::kInverted, VBB);

  // VAB.Print();
  // VBBinv.Print();

  // TMatrixD X(VAB, TMatrixD::kMult, VBBinv);
  // X *= -1.;

  TH1D *h = new TH1D("hTLS", "hTLS", fN, fTrueX1, fTrueX2);
  for (int j=0; j<fN; j++)
    h->SetBinContent(j+1, x(j));

  return h;
}

UnfoldingResult
UnfoldingUtils::UnfoldRichardsonLucy(const int nIterations)
{
  // See eq. 2.11, J. Bardsley & J. Nagy,
  // SIAM. J. MATRIX ANAL. APPL. Vol 27 No. 4, pp. 1184-1197
  // Also
  // SIAM J. SCI. COMPUT. Vol. 25, No. 4, pp. 1326–1343

  UnfoldingResult result;
  static int id = 0; id++;       // Unique ID

  result.LCurve = new TGraph(nIterations);
  result.GcvCurve = new TGraph(nIterations);
  result.LCurve->SetNameTitle(Form("LCurve_RL_%d",id),
                              Form("Richardson-Lucy L-Curve;"
                                   "Total #chi^{2};"
                                   "Total curvature of x_{k}"));
  TMatrixD A = GetA();
  TVectorD b = Getb();
  TVectorD ones = Ones(fM);

  result.XReg.ResizeTo(fN, nIterations);
  result.WReg.ResizeTo(fN, nIterations);

  // Prior vector x0 must be positive
  TVectorD x0 = Ones(fN);
  if (fHistPrior)
    x0 = Hist2Vec(fHistPrior);

  for (int j=0; j<x0.GetNrows(); j++)
  {
    if (x0(j)<=0)
    {
      Warning("UnfoldingUtils::UnfoldRichardsonLucy()",
              "Initial point x0(%d) = %g must be positive, "
              "setting to 1.0", j, x0(j));
      x0(j) = 1.;
    }
  }

  TMatrixD AT(TMatrixD::kTransposed, A);
  TVectorD AT1 = AT*ones;
  TVectorD x(x0);
  TVectorD bkg(fM);
  TVectorD bvar(fM); // Variance of meas. data
  for (int i=0; i<fM; i++)
    bvar(i) = fMatB(i,i);

  result.kLcv = 0;
  TVectorD kIter(nIterations);

  for (int k=0; k<nIterations; k++)
  {
    printf("Richardson-Lucy iteration %d\r", k+1);
    kIter(k) = k+1;

    // R-L result for iteration k
    TVectorD xd = ElemDiv(x, AT1);
    x = ElemMult(xd, AT*ElemDiv(b+bvar, A*x + bkg + bvar));
    for (int j=0; j<fN; j++)
    {
      result.WReg(j,k) = x(j);
      result.XReg(j,k) = x(j)*fVecXini(j);
    }

    // L-Curve
    TVectorD r = A*x-b;
    TVectorD chi2vec = ElemDiv(ElemMult(r,r), bvar);
    result.LCurve->SetPoint(k,chi2vec.Sum(),Curvature(x));

  } // end iteration loop
  cout << endl;

  // Assign result.Lcurvature and result.kLcv
  result.LCurvature = LogCurvature(result.LCurve, kIter, result.kLcv);
  result.LCurvature->GetXaxis()->CenterTitle();
  result.LCurvature->GetYaxis()->CenterTitle();
  result.LCurvature->GetXaxis()->SetTitleOffset(1.3);
  result.LCurvature->GetYaxis()->SetTitleOffset(1.3);
  result.lambdaLcv = result.kLcv; // Here, "lambda" is the iteration.

  result.WRegHist = Matrix2Hist(result.WReg, Form("W_RL_%d",id),
                                fTrueX1,fTrueX2,0,nIterations);
  result.WRegHist->GetXaxis()->CenterTitle();
  result.WRegHist->GetYaxis()->CenterTitle();
  result.WRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.WRegHist->GetYaxis()->SetTitleOffset(1.8);

  result.XRegHist = Matrix2Hist(result.XReg, Form("X_RL_%d",id),
                                fTrueX1,fTrueX2,0,nIterations);
  result.XRegHist->GetXaxis()->CenterTitle();
  result.XRegHist->GetYaxis()->CenterTitle();
  result.XRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.XRegHist->GetYaxis()->SetTitleOffset(1.8);

  return result;
}

UnfoldingResult
UnfoldingUtils::UnfoldPCGLS(const int nIterations,
                            int ltype,
                            TString opt,
                            const GSVDResult *gsvd,
                            const TH2 *hA,
                            const TH1 *hb,
                            const TH1 *hXini)
{
  TMatrixD L = LMatrix(fN, ltype);
  return UnfoldPCGLS(nIterations,
                     L,
                     opt,
                     gsvd,
                     hA,
                     hb,
                     hXini);
}

UnfoldingResult
UnfoldingUtils::UnfoldPCGLS(const int nIterations,
                            TMatrixD &L,
                            TString opt,
                            const GSVDResult *gsvd,
                            const TH2 *hA,
                            const TH1 *hb,
                            const TH1 *hXini)
{
  //
  // Key ingredients and their dimensions
  // -----------------------------------------------------------------
  // A: m x n  (hResp) TH2 is true (y) vs meas (x), as in RooUnfold.
  // b: m      (hb)
  // x: n      (solutions)
  // L: p x n  p has no direct restrictions
  // W: n x ?  Basis vectors of L nullspace. W != 0 for p < n.
  // hXini:    x^{ini} from Hocker et al, NIM A 372 (1996) 469-481.
  //           Optional. Improves regularization of steep functions.
  // -----------------------------------------------------------------
  //
  // Some equations below reference the book "Rank Deficient and
  // Discrete Ill-Posed Problems" by P.C. Hansen, sections 2.3.2 and
  // 6.3.
  //

  UnfoldingResult result;
  static int id = 0; id++;       // Unique ID

  TMatrixD A = GetA(opt);
  TVectorD b = Getb(opt);
  TVectorD xini(fVecXini);
  if (hA) A = Hist2Matrix(hA);
  if (hb) b = Hist2Vec(hb);
  if (hXini)
  {
    xini = Hist2Vec(hXini);
  }

  TMatrixD W = Null(L);          // Nullspace of L
  int p = L.GetNrows();

  TVectorD g2(fN);
  if (gsvd)
    g2 = ElemMult(gsvd->gamma, gsvd->gamma);
  TMatrixD B(fMatB);
  if (opt.Contains("~"))
    B.UnitMatrix();

  // Columns are filter factors at iteration k
  result.F.ResizeTo(fN,nIterations);
  TVectorD Fd(fN);
  TVectorD fk(fN); // column k of result.F
  double gcvMin = 1e99;

  result.WReg.ResizeTo(fN,nIterations);
  result.XReg.ResizeTo(fN,nIterations);
  result.WRegErr.ResizeTo(fN,nIterations);
  result.XRegErr.ResizeTo(fN,nIterations);
  result.wCov.ResizeTo(fN,fN);
  result.xCov.ResizeTo(fN,fN);

  result.LCurve = new TGraph(nIterations);
  result.GcvCurve = new TGraph(nIterations);

  result.LCurve->SetNameTitle(Form("LCurve_PCGLS_%d",id),
                              Form("PCGLS L-Curve;"
                                   "Residual norm ||Ax_{k}-b||_{2};"
                                   "Solution norm ||Lx_{k}||_{2}"));
  result.GcvCurve->SetNameTitle(Form("GcvCurve_PCGLS_%d",id),
                                Form("PCGLS GCV Curve;"
                                     "Iteration k;"
                                     "G(k)"));

  // Store q1 vectors as columns of Q1n for re-orthogonalization
  TMatrixD Q1n(p, nIterations+1);

  // For transformation to standard form
  // T = dagger(AW)*A and x0 = W*dagger(AW)*b
  TMatrixD T(0,0);
  TVectorD x0(fN); x0.Zero();
  if (W.GetNcols()>0)   // if L has a nontrivial null space
  {
    TMatrixD AW(A*W);
    TMatrixD AWinv = MoorePenroseInverse(AW);
    T.ResizeTo(AWinv.GetNrows(), fN);
    T = AWinv*A;                  // eq. 2.47
    x0 = W*AWinv*b;               // eq. 2.46
  }

  TVectorD x(x0);                    // kth solution
  TMatrixD AT(TMatrixD::kTransposed, A);
  TVectorD r = b - A*x0;
  TVectorD s = AT*r;
  TVectorD q(p), q1(p);

  LTSolve(q1, L, s);        // q1 = L11^{-T} * s. Length = p.
  LSolve(q, L, q1, W, T);   // q = L_A^+ * q1

  TVectorD z(q);
  double dq = s*q, dq2=0;
  double alpha=0, beta=0;

  // Reorthogonalize using MGS to improve accuracy
  double q1norm = TMath::Sqrt(q1*q1);
  if (q1norm > 0) q1norm = 1./q1norm;
  TMatrixDColumn(Q1n, 0) = q1norm*q1;

  for (int k=1; k<=nIterations; k++)
  {
    printf("PCGLS iteration %d\r", k);
    TVectorD Az = A*z;
    alpha = dq / Az.Norm2Sqr();
    x += alpha*z;
    r -= alpha*Az;
    s  = AT*r;
    LTSolve(q1, L, s);

    // Reorthogonalize q1 to previous q1 vectors
    for (int i=0; i<k; i++)
    {
      TVectorD qi = TMatrixDColumn(Q1n, i);
      q1 -= (qi*q1)*qi; // Modified Graham-Schmidt
      q1norm = TMath::Sqrt(q1*q1);
      if (q1norm > 0) q1norm = 1./q1norm;
      TMatrixDColumn(Q1n, k) = q1norm*q1;
    }

    LSolve(q, L, q1, W, T);
    dq2 = s*q;
    beta = dq2/dq;
    dq = dq2;
    z = q + beta*z;

    // Results
    for (int j=0; j<fN; j++)
      result.WReg(j,k-1) = x(j);
    for (int j=0; j<fN; j++)
      result.XReg(j,k-1) = x(j)*xini(j);

    // Parameter optimization analysis -------------------------------

    // L-Curve
    result.LCurve->SetPoint(k-1, TMath::Sqrt(r*r), TMath::Sqrt((L*x)*(L*x)));

    // Filter factors (if GSVD gamma values available)
    if (gsvd)
    {

      if (k==1)
      {
        fk = alpha*g2;
        Fd = g2 - ElemMult(g2, fk) + beta*g2;
        TMatrixDColumn(result.F, k-1) = fk;
      }
      else
      {
        fk += alpha*Fd;
        Fd = g2 - ElemMult(g2, fk) + beta*Fd;
        TMatrixDColumn(result.F, k-1) = fk;
      }
      if (k > 2)
      {
        for (int i=0; i<fN; i++)
        {
          if (TMath::Abs(result.F(i,k-2)-1) < 1e-4)
            result.F(i,k-1) = 1.0;
          if (TMath::Abs(result.F(i,k-3)-1) < 1e-4)
            result.F(i,k-1) = 1.0;

          if (result.F(i,k-1) > 1.)
            result.F(i,k-1) = 1.0;

          fk(i) = result.F(i,k-1);
        }
      }

      // Generalized cross-validation curve
      TVectorD resid = A*x - b;
      double gcv = TMath::Sqrt(resid*resid) / (fM - fk.Sum());
      result.GcvCurve->SetPoint(k-1, k, gcv);
      if (gcv < gcvMin)
      {
        gcvMin = gcv;
        result.kGcv = k;
      }

      // Compute A^#, the regularized inverse of A, using these filter
      // factors.
      TMatrixD FCd(fN,fN);
      TMatrixD UT(TMatrixD::kTransposed, gsvd->U);
      for (int i=0; i<fN; i++)
      {
        double a = gsvd->alpha(i);
        FCd(i,i) = (a > 0.) ? result.F(i,k-1)/a : 0.;
      }
      TMatrixD Ap = gsvd->X * FCd * UT;

      // Covariance matrices:
      // Cov(w) = Ap*B*Ap'
      // Cov(x)_ik = xini_i Cov(w)_ik xini_k
      result.wCov = Ap * TMatrixD(B, TMatrixD::kMultTranspose, Ap);
      for (int i=0; i<fN; i++)
      {
        for (int j=0; j<fN; j++)
        {
          result.xCov(i,j) = fVecXini(i) * result.wCov(i,j) * fVecXini(j);
        }
      }

      // Assign uncertainty
      for (int j=0; j<fN; j++)
        result.WRegErr(j,k-1) = result.wCov(j,j);
      for (int j=0; j<fN; j++)
        result.XRegErr(j,k-1) = result.xCov(j,j);

    }

  } // end iteration loop

  cout << endl;

  // Solutions to Atilde*w = btilde
  result.WRegHist = Matrix2Hist(result.WReg, Form("W_CG_%d",id),
                                fTrueX1, fTrueX2,0,nIterations);
  result.WRegHist->GetXaxis()->CenterTitle();
  result.WRegHist->GetYaxis()->CenterTitle();
  result.WRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.WRegHist->GetYaxis()->SetTitleOffset(1.8);
  result.WRegHist->SetTitle("CGLS solutions w_{k};w;iteration k");

  // Solutions x_j = w_j * xini_j
  result.XRegHist = Matrix2Hist(result.XReg, Form("X_CG_%d",id),
                                fTrueX1, fTrueX2,0,nIterations);
  result.XRegHist->GetXaxis()->CenterTitle();
  result.XRegHist->GetYaxis()->CenterTitle();
  result.XRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.XRegHist->GetYaxis()->SetTitleOffset(1.8);

  result.XRegHist->SetTitle("CGLS solutions x_{k};x;iteration k");

  // Add uncertainties
  TH2D *errw = Matrix2Hist(result.WRegErr, Form("errW_CG_%d",id),
                           fTrueX1, fTrueX2,0,nIterations);
  TH2D *errx = Matrix2Hist(result.XRegErr, Form("errX_CG_%d",id),
                           fTrueX1, fTrueX2,0,nIterations);

  for (int j=1; j<=fN; j++)
  {
    for (int k=1; k<nIterations; k++)
    {
      result.WRegHist->SetBinError(j,k,errw->GetBinContent(j,k));
      result.XRegHist->SetBinError(j,k,errx->GetBinContent(j,k));
    }
  }

  result.hGcv =
    result.XRegHist->ProjectionX(Form("cg_%d_iter%d",id,result.kGcv),
                                 result.kGcv,result.kGcv);

  return result;
}

TGraph *
UnfoldingUtils::LogCurvature(TGraph *g, const TVectorD &tVec, int &kMax)
{
  // Compute signed log curvature of the parametric curve
  // g(t) = (a(t), b(t))
  // where (x(t), y(t)) = (log(a(t)), log(b(t))).
  // Curvature = (x'y'' - y'x'') / (x'^2 + y'^2)^1.5
  // kMax is the index of the maximum point.

  int n = g->GetN();
  TGraph *gk = new TGraph(n-2);

  TVectorD x(n);
  TVectorD y(n);
  for (int i=0; i<n; i++)
  {
    x(i) = TMath::Log10(g->GetX()[i]);
    y(i) = TMath::Log10(g->GetY()[i]);
  }

  // Differentiation operators
  // TMatrixD D1 = LMatrix(n,k1DerivNoBC);
  // TMatrixD D2 = LMatrix(n,k2DerivNoBC);

  // // Approximate derivatives
  // TVectorD x1 = D1*x;
  // TVectorD x2 = D2*x;
  // TVectorD y1 = D1*y;
  // TVectorD y2 = D2*y;

  TVectorD x1(n-1);
  TVectorD x2(n-2);
  TVectorD y1(n-1);
  TVectorD y2(n-2);

  // First derivatives
  for (int i=0; i<n-1; i++)
  {
    x1(i) = (x(i+1) - x(i)) / (tVec(i+1) - tVec(i));
    y1(i) = (y(i+1) - y(i)) / (tVec(i+1) - tVec(i));
  }
  // Second derivatives
  for (int i=0; i<n-2; i++)
  {
    x2(i) = (x1(i+1) - x1(i)) / (tVec(i+1) - tVec(i));
    y2(i) = (y1(i+1) - y1(i)) / (tVec(i+1) - tVec(i));
  }

  // Compute curvature vs. t
  double cmax = -1e99;
  for (int i=0; i<n-2; i++)
  {
    double numer = x1(i)*y2(i) - x2(i)*y1(i);
    double denom = TMath::Power(x1(i)*x1(i) + y1(i)*y1(i), 1.5);
    double c = numer / denom;
    gk->SetPoint(i, tVec(i), c);

    // Find maximum curvature and its array position
    if (c>cmax)
    {
      cmax = c;
      kMax = i;
    }

  }
  return gk;
}

TH1D *
UnfoldingUtils::XHist(TVectorD &x, TString base, int k,
                      double xMin, double xMax,
                      double normto, TString /*opt*/)
{
  // base should have a form like "CGLS"+"id"
  // k is the iteration step
  // The factor normto could be hb->Integral();
  const char *name  = Form("h%s_%d", base.Data(), k);
  const char *title = Form("%s (%d iterations)", base.Data(), k);
  TH1D *hx = Vec2Hist(x, xMin, xMax, name, title);

  if (normto != 0.0)
  {
    double norm = hx->Integral() ? normto/hx->Integral() : 0;
    hx->Scale(norm);
  }
  // TODO assign proper stat uncertainty
  return hx;
}

TMatrixD
UnfoldingUtils::LMatrix(const int n, const int kind, double eps)
{
  // Return p x n smoothing matrix. The number of rows p is
  // case-dependent.

  // Top row and left column.
  // Make sure that row[0] = col[0]
  double row[n];
  for (int j=0; j<n; j++) row[j]=0.0;

  if (kind == kUnitMatrix)   // n x n identity matrix. p = n. W=0.
  {
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    return Toeplitz(n,n,col,row);
  }
  if (kind == k1DerivNoBC)   // 1st deriv, no BC assumptions. p = n-1. W=const.
  {
    double col[n-1];
    for (int i=0; i<n-1; i++) col[i]=0.0;
    row[0] = col[0] = -1;
    row[1] = 1;
    return Toeplitz(n-1,n,col,row);
  }
  if (kind == k2DerivNoBC)   // 2nd deriv, no BC assumptions. p = n-2. W=const, linear.
  {
    double col[n-2];
    for (int i=0; i<n-2; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    row[1] = -2;
    row[2] = 1;
    return Toeplitz(n-2,n,col,row);
  }
  if (kind == k1DerivBC0)   // 1st deriv, BC=0 on L and R. p = n+1. W=0.
  {
    double col[n+1];
    for (int i=0; i<n+1; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    col[1] = -1;
    return Toeplitz(n+1,n,col,row);
  }
  if (kind == k2DerivBC0)   // 2nd deriv, BC=0 on L and R. p = n. W=0.
  {
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = -2+eps;
    row[1] = col[1] = 1;
    return Toeplitz(n,n,col,row);
  }
  if (kind == k1DerivBCR)   // 1st deriv, reflect at L,R. p = n-1. W=const. Same as k1DerivNoBC
  {
    double col[n-1];
    for (int i=0; i<n-1; i++) col[i]=0.0;
    row[0] = col[0] = -1;
    row[1] = 1;
    return Toeplitz(n-1,n,col,row);
  }
  if (kind == k2DerivBCR)   // 2nd deriv, reflect at L,R. p = n. W=const.
  {
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = -2+eps;
    row[1] = col[1] = 1;
    TMatrixD L = Toeplitz(n,n,col,row);
    L(0,0) = -1+eps;
    L(n-1,n-1) = -1+eps;
    return L;
  }

  TMatrixD Ldefault(n,n);
  Ldefault.UnitMatrix();
  return Ldefault;
}

void
UnfoldingUtils::SetTH1Props(TH1 *h,
                            Int_t linecolor,
                            Int_t fillcolor,
                            Int_t markercolor,
                            Int_t markerstyle,
                            Double_t markersize)
{
  h->SetLineColor(linecolor);
  h->SetFillColor(fillcolor);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  h->SetMarkerSize(markersize);
}

TGraph *
UnfoldingUtils::ResidualNorm(TObjArray *hists, double stepSize)
{
  static int id=0; id++;

  // Check for Ax-b histos in hists
  TObjArray *subList = new TObjArray();
  for (int i=0; i<hists->GetEntries(); i++)
  {
    TObject *obj = hists->At(i);
    TString name = obj->GetName();
    TString cl   = obj->ClassName();
    if (name.Contains("hr") && cl.Contains("TH1"))
    {
      subList->Add(obj);
      if (0)
        Info("UnfoldingUtils::ResidualNorm()", "Added %s",name.Data());
    }
  }
  int np = subList->GetEntries();
  TGraph *g = new TGraph(np);
  for (int i=0; i<np; i++)
  {
    TH1 *h = (TH1 *)subList->At(i);
    TVectorD r = Hist2Vec(h);
    g->SetPoint(i, stepSize*(i), TMath::Sqrt(r*r));
  }

  if (g->GetN()==0)
    Warning("UnfoldingUtils::ResidualNorm()", "No points in graph. np=%d", np);

  g->SetLineColor(kBlue-2);
  g->SetMarkerColor(kBlue-2);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1.);
  g->SetLineWidth(2);
  g->SetNameTitle(Form("rn_%d",id),
                  Form("Residual Norm;"
                       "#lambda;"
                       "Residual norm ||Ax^{(#lambda)}-b||_{2};"));
  return g;
}

TH2D *
UnfoldingUtils::UnfoldCovMatrix(int nTrials,
                                int algo,
                                double regPar,
                                TString opt,
                                TObjArray *bHists)
{
  // Propagate the measured data cov matrix to the unfolded result.
  static int id=0; id++;

  if (bHists && nTrials > bHists->GetEntries())
  {
    Error("UnfoldingUtils::UnfoldCovMatrix()",
          "nTrials (%d) > number of trial data sets (%d)",
          nTrials, bHists->GetEntries());
    return 0;
  }

  // seed = 0 signals to use the TUUID identifier
  int seed = 0;
  TRandom3 rand(seed);

  // Unfolded result
  TH1D *hUnf = 0;

  // For iterative methods, regPar serves as the number of iterations.
  // int nk = regPar;

  // Measurement covariance
  TMatrixD B(fMatB);
  if (opt.Contains("~"))
    B = LMatrix(fM, kUnitMatrix); // B~ is unit matrix

  // TDecompChol decomp(fMatB);
  // TMatrixD UT(TMatrixD::kTransposed, decomp.GetU());

  // Returned covariance matrix
  TString title(Form("Covariance matrix for %s result",
                     Algorithm(algo).Data()));
  TH2D *hcov = new TH2D(Form("hcov%d",id),title.Data(),
                        fN,fTrueX1,fTrueX2,fN,fTrueX1,fTrueX2);

  // Gaussian perturbation
  TVectorD g(fM);

  // Measured data + MC perturbation for one pseudo-experiment
  TH1D *hbTrial = (TH1D *)fHistMeas->Clone(Form("hbTrial%d",id));

  // Mean result of unfolding MC trials E[x]
  TH1D *hxMean = new TH1D(Form("hxMean%d",id),"mean",fN,fTrueX1,fTrueX2);

  // Loop 1: Compute E[x] as hxMean
  for (int t=0; t<nTrials; t++)
  {

    cout << Form("UnfoldCovMatrix(): computing mean result %d\r",t) << flush;

    if (bHists)
    {
      fHistMeas = (TH1D *)bHists->At(t);
      fVecb = Hist2Vec(fHistMeas);
      fVecbErr = Hist2Vec(fHistMeas, "unc");
      fMatATilde = DivColsByVector(fMatA, fVecbErr);
      fVecbTilde = ElemDiv(fVecb, fVecbErr);
    }

    else
    {
      for (int i=0; i<fM; i++)
      {
        g(i) = rand.Gaus(0,1);
      }

      g *= B;//UT;
      for (int i=0; i<fM; i++)
      {
        hbTrial->SetBinContent(i+1, fVecb(i) + g(i));
        hbTrial->SetBinError(i+1, fHistMeas->GetBinError(i+1));
        fVecb(i) = hbTrial->GetBinContent(i+1);
        fVecbErr(i) = hbTrial->GetBinError(i+1);
      }
    }

    switch (algo)
    {
    case kGSVDAlgo:
      hUnf = GSVDAnalysis(fMatL, regPar, 0, 0, opt)->xregHist;
      break;
    case kSVDAlgo:
      hUnf = UnfoldSVD(regPar, 0, opt, 0, 0, 0);
      break;

      // case kRichLucyAlgo:
      //   hUnf = UnfoldRichardsonLucy(nk);
      //   break;

      // case kChi2MinAlgo:
      //  h = UnfoldChi2Min();
      //  break;
      // TODO: Implement hLbest or hGcv for PCGLS
      // case kPCGLSAlgo:
      //   hUnf = UnfoldPCGLS(nk,fMatL,opt,0,0,hbTrial);
      //   break;
    }

    for (int j=1; j<=fN; j++)
      hxMean->AddBinContent(j,hUnf->GetBinContent(j)/nTrials);
  }
  cout << endl;
  rand.SetSeed(seed);

  // Loop 2: Compute E[(x-E[x])*(x-E[x])'] as hcov
  for (int t=0; t<nTrials; t++)
  {

    cout << Form("UnfoldCovMatrix(): computing covariance %d\r",t) << flush;

    if (bHists)
    {
      fHistMeas = (TH1D *)bHists->At(t);
      fVecb = Hist2Vec(fHistMeas);
      fVecbErr = Hist2Vec(fHistMeas, "unc");
      fMatATilde = DivColsByVector(fMatA, fVecbErr);
      fVecbTilde = ElemDiv(fVecb, fVecbErr);
    }
    else
    {
      for (int i=0; i<fM; i++)
      {
        g(i) = rand.Gaus(0,1);
      }
      g *= B;//UT;
      for (int i=0; i<fM; i++)
      {
        hbTrial->SetBinContent(i+1, fVecb(i) + g(i));
        hbTrial->SetBinError(i+1, fHistMeas->GetBinError(i+1));
      }
    }
    switch (algo)
    {
    case kGSVDAlgo:
      hUnf = GSVDAnalysis(fMatL, regPar, 0, 0, opt)->xregHist;
      break;
    case kSVDAlgo:
      hUnf = UnfoldSVD(regPar, 0, opt, 0, 0, 0);
      break;
    }

    double dxj=0, dxk=0, cjk=0;
    for (int j=1; j<=fN; j++)
    {
      for (int k=1; k<=fN; k++)
      {
        dxj = hUnf->GetBinContent(j) - hxMean->GetBinContent(j);
        dxk = hUnf->GetBinContent(k) - hxMean->GetBinContent(k);
        cjk = dxj*dxk/(nTrials-1);
        cjk += hcov->GetBinContent(j,k);
        hcov->SetBinContent(j,k,cjk);
      }
    }

  }
  cout << endl;

  // Should ComputeRescaledSystem() be called here to reset back to original inputs?
  return hcov;
}

// -------------------------------------------------------------------
// ----------------------- Support classes ---------------------------
// -------------------------------------------------------------------

GSVDResult::GSVDResult() :
  n(0),                 // Column count of A (or L)
  m(0),                 // Row count of A
  p(0),                 // Row count of L
  lambda(0.0),         // Regularization parameter used
  UHist(0),           // Left sing. vectors
  XHist(0),          // GSVD basis vectors
  wregHist(0),        // scaled result w^lambda
  xregHist(0),        // xini_j * w^lambda_j (solution)
  bregHist(0),        // Refolded distribution A*xreg
  UTbAbs(0),          // Vector of |U'_i*b| values
  coeffAbs(0),        // Vector of |U'_i*b|/alpha_i values
  regcAbs(0)         // Regularized (filtered) coeffs.
{
}
