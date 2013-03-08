#include "UnfoldingUtils.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TDecompSVD.h"
#include "TDecompQRH.h"
#include "TDecompChol.h"
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
#endif

using std::cout;
using std::endl;

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

UnfoldingUtils::UnfoldingUtils(TH2D* hA, 
			       TH1D* hMeas, 
			       TH2D* hMeasCov, 
			       TH1D* hXini, 
			       TH1D* hXtrue, 
			       TH1D* hEff) :
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

  if (!BinningOk()) {
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

  if (fHistXini) {
    fVecXini = Hist2Vec(fHistXini);
  }
  else {
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
  if (!fHistAProb) {
    fHistAProb = (TH2D*) fHistA->Clone("fHistAProb");
    NormalizeXSum(fHistAProb, fHistEff);
    fHistAProb->SetTitle("Probability matrix #hat{A}");
  }
  fMatAhat = Hist2Matrix(fHistAProb);
  
  // Data uncertainty
  for (int i=0; i<fM; i++)
    fVecbErr(i) = fHistMeas->GetBinError(i+1);
  
  // Create rescaled (~) quantities
  if (fHistMeasCov) {
    
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
  else { 
    // The usual case - b has indep. errors
    for (int i=0; i<fM; i++) {
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

  for (int i=0;i<fM;i++)
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
  if (nMeas != fM) {
    gROOT->Warning("UnfoldingUtils::BinningOk()", 
		   "Meas. bin mismatch: hMeas %d, TH2 (x-axis) %d",
		   nMeas, fM);
    isok = false;
  }
  if (fHistXini) {
    int nXini = fHistXini->GetNbinsX(); // Must equal fN
    if (nXini != fN) {
      gROOT->Warning("UnfoldingUtils::BinningOk()", 
		     "True bin mismatch: x^ini %d, TH2 (y-axis) %d",
		     nXini, fN);
      isok = false;
    }
  }
  if (fHistXtrue) {
    int nXtrue = fHistXtrue->GetNbinsX(); // Must equal fN
    if (nXtrue != fN) {
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
  switch (type) {
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
UnfoldingUtils::SetPrior(TH1D* h)
{
  fHistPrior = (TH1D*)h->Clone("fHistPrior");
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

TH1D*
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
UnfoldingUtils::SVDAnalysis(TH2* hA, TH1* hb, TString opt)
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
  if (hA) {
    A.ResizeTo(hA->GetNbinsX(), hA->GetNbinsY());
    A = Hist2Matrix(hA);
  }
  if (hb) {
    b.ResizeTo(hb->GetNbinsX());
    b = Hist2Vec(hb);
  }
  if (b.GetNrows()==0) {
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

GSVDResult* 
UnfoldingUtils::GSVDAnalysis(TMatrixD& L, double lambda, TH2* hA, TH1* hb, TString opt)
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
  GSVDResult* gsvd = new GSVDResult();

  // Use stored members: 
  // Select A, \hat{A}, or \tilde{A} ("", "^", or "~")
  TMatrixD A = GetA(opt);
  TVectorD b = Getb(opt);

  // Or, use passed-in histograms (supercedes opt)
  if (hA) {
    A.ResizeTo(hA->GetNbinsX(), hA->GetNbinsY());
    A = Hist2Matrix(hA);
  }
  if (hb) {
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
  for (int i=0; i<n; i++) {

    if (i >= n-p) {
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
  for (int i=0; i<n; i++) {
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
  for (int i=0; i<n; i++) {
    for (int k=0; k<n; k++) {
      gsvd->covx(i,k) = fVecXini(i) * gsvd->covw(i,k) * fVecXini(k);
    }
  }

  // Compute cov(x)^-1 (without using cov(x), see Hocker eq. 53)
  double AAjk = 0;
  for (int j=0; j<n; j++) {
    for (int k=0; k<n; k++) {
      AAjk = 0;
      for (int i=0; i<n; i++) {
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
  for (int j=0; j<n; j++) {
    gsvd->wregHist->SetBinError(j+1, TMath::Sqrt(gsvd->covw(j,j)));
    gsvd->xregHist->SetBinError(j+1, TMath::Sqrt(gsvd->covx(j,j)));
  }
  
  // Absolute values for plotting
  TVectorD utbAbs(utb);
  TVectorD cAbs(c);
  TVectorD rcAbs(ElemMult(f,utb));   // TVectorD rcAbs(regc);
  for (int i=0; i<n; i++) {
    if (utbAbs(i) < 0) utbAbs(i) *= -1;
    if (cAbs(i) < 0)     cAbs(i) *= -1;
    if (rcAbs(i) < 0)   rcAbs(i) *= -1;
  }
  gsvd->UTbAbs   = Vec2Hist(utbAbs, 0,n,Form("gsvd_utb_abs%d",id), "#||{u^{T}#upointb} ");
  gsvd->coeffAbs = Vec2Hist(cAbs, 0,n,Form("gsvd_c_abs%d",id), "#||{u^{T}#upointb}/#alpha ");
  gsvd->regcAbs  = Vec2Hist(rcAbs, 0,n,Form("gsvd_rc_abs%d",id), "f#||{u^{T}#upointb} ");

  SetTH1Props(gsvd->UTbAbs,   kBlue, 0, kBlue, kFullSquare, 1.0);
  SetTH1Props(gsvd->coeffAbs, kRed, 0, kRed, kOpenSquare, 1.0);
  SetTH1Props(gsvd->regcAbs, kMagenta+1, 0, kMagenta+1, kOpenCircle, 1.0);
  SetTH1Props(gsvd->xregHist, kGreen+2, 0, kGreen+2, kFullCircle, 1.0);

  return gsvd;
}

TCanvas* 
UnfoldingUtils::DrawSVDPlot(SVDResult svdhists, double ymin, double ymax, TString opt)
{
  static int i=0; i++;
  TCanvas* c = new TCanvas(Form("csvd%d",i), Form("csvd%d",i), 1);

  // Draw frame histogram to set limits, title, etc.
  int nx = svdhists.sigma->GetNbinsX();
  TH1F* hsvd = new TH1F(Form("hsvd%d",i), "SV Components;column index i;", 
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
  
  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
  if (opt.Contains("sig"))
    leg->AddEntry(svdhists.sigma, svdhists.sigma->GetTitle(), "p");
  leg->AddEntry(svdhists.UTb, svdhists.UTb->GetTitle(), "ep");
  leg->AddEntry(svdhists.coeff, svdhists.coeff->GetTitle(), "ep");
  leg->SetFillColor(kNone);
  leg->Draw();
  
  return c;
}

TCanvas* 
UnfoldingUtils::DrawGSVDPlot(GSVDResult* gsvd, double ymin, double ymax, TString opt)
{
  static int i=0; i++;
  TCanvas* c = new TCanvas(Form("cgsvd%d",i), Form("cgsvd%d",i), 1);
  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99,
			     Form("NDF_{eff} = %.1f", 
				  gsvd->f.Sum()));
  int nx = gsvd->UTbAbs->GetNbinsX();
  TH1F* h = new TH1F(Form("hgsvd%d",i), 
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
  if (opt.Contains("c")) {
    gsvd->coeffAbs->Draw("plsame");
    leg->AddEntry(gsvd->coeffAbs, gsvd->coeffAbs->GetTitle(), "ep");
  // TODO: add damped coeffs f|U'*b|/alpha
  }
  
  leg->SetFillColor(kNone);
  leg->Draw();

  if (opt.Contains("~")) {
    TLine one;
    TAxis* ax = h->GetXaxis();
    one.DrawLine(ax->GetXmin(), 1.0, ax->GetXmax(), 1.0);
  }

  return c;
}

TestProblem
UnfoldingUtils::MonteCarloConvolution(const int m, 
				      const int n,
				      const double xm1,
				      const double xm2,
				      const double xt1,
				      const double xt2,
				      TF1* truthFn, 
				      TF1* kernelFn, 
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
  t.Response = Matrix2Hist(R, Form("R%d",id), xm1, xm2, xt1, xt2);
  t.Response->SetTitle(Form("%d x %d convolution matrix;"
			    "s (observed);t (true)", m, n));

  TH2D* RMC = new TH2D(Form("R_MC%d",id), "A_{MC}", 
			m,xm1,xm2,n,xt1,xt2);

  // Model a true and a measured distribution
  // There is no bIdeal for this problem.
  t.xTruth    = new TH1D("hTrue",    "",         n, xt1, xt2);  
  t.xTruthEst = new TH1D("hTrueEst", "hTrueEst", n, xt1, xt2);  
  t.bNoisy    = new TH1D("hMeas",    "hMeas",    m, xm1, xm2);  
  t.xTruthEst->Sumw2();
  t.bNoisy->Sumw2();
  for (Int_t i=0; i<nEvents; i++) {
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
  for (int j=1; j<=n; j++) {
    double val = truthFn->Eval(t.xTruth->GetBinCenter(j));
    t.xTruth->SetBinContent(j, val);
  }
  t.xTruth->Scale(nEvents/t.xTruth->Integral());

  return t;
}

TestProblem
UnfoldingUtils::ShawSystem(const int n, double noise)
{
  TMatrixD A(n,n);
  TVectorD x(n);
  TVectorD b_ideal(n);
  ShawSystem(n, A, x, b_ideal, 0.);
  TVectorD b(b_ideal);

  // Add Gaussian white noise to b
  if (noise > 0.) {
    TRandom3 r3;
    for (int j=0; j<n; j++) b(j) += noise*r3.Gaus();
  }
  // Assign output struct members.
  // No xTruthEst histogram for this problem. Don't ask for it!
  TestProblem t;
  t.Response = Matrix2Hist(A, "Shaw_A",0.,1.,0.,1.);
  t.xTruth = Vec2Hist(x, 0., 1., "Shaw_x","Truth x fn.");
  t.bIdeal = Vec2Hist(b_ideal, 0., 1., "Shaw_b_ideal","Meas. b fn.");
  t.bNoisy = Vec2Hist(b, 0., 1., "Shaw_b","Meas. b fn.");
  for (int j=0; j<n; j++)
    t.bNoisy->SetBinError(j+1, noise);
  return t;
}

void 
UnfoldingUtils::ShawSystem(const int n, TMatrixD& A, TVectorD& x, TVectorD& b, double noise)
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
  
  if (n%2) {
    gROOT->Error("UnfoldingUtils::ShawSystem()", "Even binning required");
    return;
  }

  // For computing A
  double h = TMath::Pi()/n;
  double co[n], psi[n];

  // For computing x vector
  double a1=2, a2=1, c1=6, c2=2, t1=0.8, t2=-0.5;

  for (int i=0; i<n; i++) {
    co[i] = TMath::Cos(-TMath::PiOver2() + (0.5 + i)*h);
    psi[i] = TMath::Pi()*TMath::Sin(-TMath::PiOver2() + (0.5 + i)*h);

    double arg1 = -c1*TMath::Power(-TMath::PiOver2() + (0.5 + i)*h - t1, 2);
    double arg2 = -c2*TMath::Power(-TMath::PiOver2() + (0.5 + i)*h - t2, 2);
    x(i) = a1*TMath::Exp(arg1) + a2*TMath::Exp(arg2);
  }

  for (int i=0; i<n/2; i++) {
    for (int j=i; j<n-i; j++) {
      double ss = psi[i] + psi[j];
      A(i,j) = TMath::Power((co[i]+co[j])*TMath::Sin(ss)/ss, 2);
      A(n-j-1,n-i-1) = A(i,j);
    }
    A(i,n-i-1) = TMath::Power(2*co[i], 2);
  }

  // Create upper triangular matrix from A
  TMatrixD Atri(A);
  for (int i=0; i<n; i++) {
    for (int j=i; j<n; j++) {
      if (j<=i)
  	Atri(i,j)=0.;
    }
  }
  Atri.T();

  A += Atri;
  A *= h;

  b = A*x;

  // Add Gaussian white noise to b
  if (noise > 0.) {
    TRandom3 r3;
    for (int j=0; j<n; j++) b(j) += noise*r3.Gaus();
  }
  
  return;
}

TMatrixD 
UnfoldingUtils::Hist2Matrix(const TH2* h)
{
  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();
  TMatrixD m(nx, ny);
  for (Int_t j=0; j<ny; j++) {
    for (Int_t i=0; i<nx; i++) {
      m(i,j) = h->GetBinContent(i+1,j+1);
    }
  }
  return m;
}

TH2D* 
UnfoldingUtils::Matrix2Hist(TMatrixD& A, TString hName, 
			    double x1, double x2, double y1, double y2)
{
  int m = A.GetNrows();
  int n = A.GetNcols();
  TH2D* h = new TH2D(hName.Data(),hName.Data(),m,x1,x2,n,y1,y2);
  
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      h->SetBinContent(i+1, j+1, A(i,j));
    }
  }
  
  return h;
}

TH2D* 
UnfoldingUtils::Matrix2Hist(TMatrixD& A, TString hName, 
			    double xbins[], double ybins[])
{
  int m = A.GetNrows();
  int n = A.GetNcols();
  TH2D* h = new TH2D(hName.Data(),hName.Data(),m,xbins,n,ybins);
  
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
	h->SetBinContent(i+1, j+1, A(i,j));
    }
  }
  
  return h;
}

TH2D* 
UnfoldingUtils::Matrix2Hist(TMatrixD& A, TMatrixD& errA, TString hName, 
			    double xbins[], double ybins[])
{
  int m = A.GetNrows();
  int n = A.GetNcols();
  TH2D* h = new TH2D(hName.Data(),hName.Data(),m,xbins,n,ybins);
  
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      h->SetBinContent(i+1, j+1, A(i,j));
      h->SetBinError(i+1, j+1, errA(i,j));
    }
  }
  
  return h;
}


TH1D* 
UnfoldingUtils::Vec2Hist(const TVectorD& v, Double_t x1, Double_t x2, 
			 TString name, TString title)
{
  int nb = v.GetNoElements();
  TH1D* h = new TH1D(name.Data(), title.Data(), nb, x1, x2);

  for (int i=0; i<nb; i++) {
    h->SetBinContent(i+1, v(i));
  }
  return h;
}

TVectorD 
UnfoldingUtils::Hist2Vec(const TH1* h, TString opt)
{
  // Returns TVectorD of the bin contents of the input histogram
  int nb = h->GetNbinsX();
  TVectorD v(nb);
  if (!h) return v;
  double val = 0.;
  for (Int_t i= 0; i<nb; i++) {
    if (opt.Contains("unc"))
      val = h->GetBinError(i+1);
    else
      val = h->GetBinContent(i+1);
    v(i) = val;
  }
  return v;
}

double
UnfoldingUtils::SmoothingNorm(TVectorD& x, int regtype)
{
  double sn = 0;
  switch (regtype) {
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
UnfoldingUtils::Curvature(TVectorD& x)
{
  // Eq. (38), NIM A 372 (1996) 469-481. 
  double delta=0, val=0;
  for (int i=1; i<fN-1; i++) {
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
  if (fTilde) {
    TVectorD resid = fMatATilde*x - fVecbTilde;
    chi2 = resid*resid;
  }
  else {
    TVectorD resid = fMatAhat*x - fVecb;
    chi2 = resid*(fMatBinv*resid);
  }
  
  // Additive chi^2 modifier (reg. penalty value)  
  double mod = fRegWeight*SmoothingNorm(x, fRegType);
  
  chi2 += mod;
  return chi2;
}

UnfoldingResult
UnfoldingUtils::UnfoldChiSqMin(TVectorD& regWts, TString opt)
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
  for (int j=0; j<=fN; j++) {
    xbins[j] = fTrueX1 + j*(fTrueX2-fTrueX1)/fN;
  }
  for (int k=0; k<nRegWts; k++) {
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
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(this, &UnfoldingUtils::RegChi2, nPars);
  min->SetFunction(f);
  
  // Initialize tmx array and pass to minimizer.
  double* tmx = new double[fN];
  double stepSize = 0.1;

  for (int k=0; k<nRegWts; k++) {
    for (int j=0; j<fN; j++) {
      tmx[j] = fHistPrior->GetBinContent(j+1);
      
      // Require all parameters to have a minimum positive value
      if (tmx[j] < 0) {
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
    for (int j=0; j<nPars; j++) {
      tmx[j] = min->X()[j];
      double val = tmx[j];
      double err = val < 1 ? 1. : min->Errors()[j];
      w(j) = val;
      
      result.XRegHist->SetBinContent(j+1,k+1,val);
      result.XRegHist->SetBinError(  j+1,k+1,err);

    }
    fRegWeight = 0.;
    result.LCurve->SetPoint(k, RegChi2(tmx), Curvature(w));    
  } // k loop
  
  return result;
}

TMatrixD
UnfoldingUtils::RegularizedInverseResponse(GSVDResult* gsvd, double lambda)
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
  for (int i=0; i<n; i++) {
    if (i >= n-p) {
      double g2 = gsvd->gamma(i)*gsvd->gamma(i); 
      f(i) = g2 / (g2 + lambda*lambda);
    }
    else
      f(i) = 1.0;
  }
  for (int i=0; i<fN; i++) {
    FCd(i,i) = (gsvd->alpha(i) > 0.) ? f(i)/gsvd->alpha(i) : 0.;
  }
  Ap = gsvd->X * FCd * UT; 
  return Ap;
}

UnfoldingResult
UnfoldingUtils::UnfoldTikhonovGSVD(GSVDResult* gsvd,
				   TVectorD& lambda, 
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

  // Compute the number of significant GSVD coefficients
  int nSignificantGSVDCoeffs = 0;
  double errThreshold = 2*GetbErrMean();
  for (int i=0; i<fN; i++) {
    TH1D* h = gsvd->UTbAbs;
    if (h->GetBinContent(i+1) > errThreshold)
      nSignificantGSVDCoeffs++;
    else
      break;
  }
  Printf("# coeffs > %f = %d", errThreshold, nSignificantGSVDCoeffs);

  // Stuff for computing covariance
  TMatrixD B(gsvd->covb); // Error matrix of b
  TMatrixD FCd(n,n);      // Filter factor matrix * pseudoinverse(diag(alpha))
  // TMatrixD Ap(n,m);       // Regularized inverse A^#
  TMatrixD UT(TMatrixD::kTransposed, gsvd->U);  

  // Scan over lambda values, generate nk solutions
  for (int k=0; k<nk; k++) {
    
    // Create Tikhonov filter factors for this lambda
    double l = lambda(k);
    TVectorD f(n);
    for (int i=0; i<n; i++) {
      if (i >= n-p) {
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
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
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
    for (int j=0; j<n; j++) {
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
    
    for (int j=0; j<n; j++) {
      double vxreg = covx(j,j);
      double vxinv = covxInv(j,j);
      //      double vxinv = (gsvd->covxInv)(j,j);
      double prod = vxreg*vxinv;

      // if (prod == 0 || prod < 1.) 
      // 	rho_j = 0;
      // else
      // 	rho_j = TMath::Sqrt(1 - 1./prod);

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

    if (fsum > nSignificantGSVDCoeffs) {
      result.lambdaStf = lambda(k);
      result.kStf = k;
    }
    if (gcv < gcvMin) {
      gcvMin = gcv;
      result.lambdaGcv = lambda(k);
      result.kGcv = k;
    }
    if (rho < rhoMin) {
      rhoMin = rho;
      result.lambdaRho = lambda(k);
      result.kRho = k;
    }

  }

  // Assign result.Lcurvature and result.kLcv
  result.LCurvature = LogCurvature(result.LCurve, lambda, result.kLcv);
  result.lambdaLcv = lambda(result.kLcv);

  int kg = result.kGcv+1;
  result.RhoCurve->SetName("gsvd_rho");
  result.RhoCurve->SetTitle("GSVD mean global correlation coefficients;"
			    "#lambda;#LT#rho#GT");
  result.LCurve->SetTitle("GSVD L-Curve;||Ax_{#lambda}-b||_{2};||Lx_{#lambda}||_{2}");

  result.LCurvature->SetName("gsvd_lcc");
  result.LCurvature->SetTitle("GSVD L-Curve log curvature;#lambda;log curvature");
  result.GcvCurve->SetNameTitle("gsvd_gcv","GSVD cross-validation curve;"
				"#lambda;G(#lambda)");
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

  double kbins[nk+1], xbins[fN+1];
  for (int j=0; j<=fN; j++) {
    xbins[j] = fTrueX1 + j*(fTrueX2-fTrueX1)/fN;
  }
  for (int k=0; k<nk; k++) {
    kbins[k] = lambda(k);
  }
  kbins[nk] = 2*kbins[nk-1]-kbins[nk-2];

  result.WRegHist = Matrix2Hist(result.WReg, WRegErr, Form("wregHist_%d",id), xbins, kbins);
  result.WRegHist->GetXaxis()->CenterTitle();
  result.WRegHist->GetYaxis()->CenterTitle();
  result.WRegHist->GetXaxis()->SetTitleOffset(1.8);
  result.WRegHist->GetYaxis()->SetTitleOffset(1.8);
  result.WRegHist->SetTitle("GSVD solutions: w_{#lambda};w;#lambda");


  result.XRegHist = Matrix2Hist(result.XReg, XRegErr, Form("xregHist_%d",id), xbins, kbins);
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

TH1D* 
UnfoldingUtils::UnfoldSVD(double lambda, 
			  TObjArray* output,
			  TString opt,
			  TH2* hA,
			  TH1* hb,
			  TH1* hXini)
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

  if (fTilde) {
    A = fMatATilde;
    b = fVecbTilde;
  }
  else {
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
  for (int i=0; i<nd; i++) {
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
  for (int i=0; i<nd; i++) {
    for (int k=0; k<nd; k++) {
      Xcov(i,k) = xini(i)*Wcov(i,k)*xini(k);
    }
  }

  // X inverse (without using Xcov, see eq. 53)
  TMatrixD Xinv(nd, nd);
  double AAjk = 0;
  for (int j=0; j<nd; j++) {
    for (int k=0; k<nd; k++) {
      AAjk = 0;
      for (int i=0; i<nd; i++) {
	AAjk += A(i,j)*A(i,k);
      }
      Xinv(j,k) = AAjk/xini(j)/xini(k);
    }
  }
  
  // Add components to output list
  if (output) {
    TH1D* hs = Vec2Hist(s_vec, 0., nd, Form("hs%d",id), "s_{i} ");
    TH1D* hd = Vec2Hist(absd,  0., nd, Form("hd%d",id), Form("#||{d_{i}} = #||{(U^{T}b)_{i}}" ));
    TH1D* hl = Vec2Hist(dlam,  0., nd, Form("hl%d",id), Form("#||{d^{(#lambda)}_{i}}, #lambda = %g ",lambda));
    TH1D* hw = Vec2Hist(w,     0., nd, Form("hw%d",id), Form("w^{#lambda = %g} ",lambda));
    TH1D* ht = Vec2Hist(tf,    0., nd, Form("ht%d",id), Form("s_{i}^{2}/(s_{i}^{2}+#lambda^{2}), #lambda = %g ",lambda));
    TH1D* hr = Vec2Hist(resid, 0., nd, Form("hr%d",id), Form("Ax^{(#lambda)}-b, #lambda = %g ",lambda));

    SetTH1Props(hs, kBlack, 0, kBlack, kFullCircle, 1.0);
    SetTH1Props(hd, kBlue,  0, kBlue,  kFullSquare, 1.0);
    SetTH1Props(hl, kMagenta+1,  0, kMagenta+1,  kOpenSquare, 1.0);
    SetTH1Props(hw, kCyan+2,0, kCyan+2,kOpenCircle, 1.0);
    SetTH1Props(ht, kBlack, 0, kBlack, kOpenCircle, 1.0);

    // Covariance matrix of solution and its inverse
    TH2D* hWcov = Matrix2Hist(Wcov, Form("hWcov%d",id), 0, nd, 0, nd);
    TH2D* hXcov = Matrix2Hist(Xcov, Form("hXcov%d",id), 0, nd, 0, nd);
    TH2D* hXinv = Matrix2Hist(Xinv, Form("hXinv%d",id), 0, nd, 0, nd);

    for (int i=0; i<nd; i++) {
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
  
  TH1D* hx = Vec2Hist(x,fTrueX1,fTrueX2,Form("hSVD%d",id),
		      Form("#lambda = %g", lambda));

  for (int i=0; i<nd; i++) {
    hx->SetBinError(i+1, TMath::Sqrt(Xcov(i,i)));
  }
  
  return hx;
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

  for (int j=0; j<x0.GetNrows(); j++) {
    if (x0(j)<=0) {
      Warning("UnfoldingUtils::UnfoldRichardsonLucy()",
	      "Initial point x0(%d) = %g must be positive, "
	      "setting to 1.0", j, x0(j) );
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
  
  for (int k=0; k<nIterations; k++) {
    printf("Richardson-Lucy iteration %d\r", k+1);
    
    // R-L result for iteration k
    TVectorD xd = ElemDiv(x, AT1);
    x = ElemMult(xd, AT*ElemDiv(b+bvar, A*x + bkg + bvar));
    for (int j=0; j<fN; j++) {
      result.WReg(j,k) = x(j);
      result.XReg(j,k) = x(j)*fVecXini(j);
    }
    
    // L-Curve    
    TVectorD r = A*x-b;
    TVectorD chi2vec = ElemDiv(ElemMult(r,r), bvar);
    result.LCurve->SetPoint(k,chi2vec.Sum(),Curvature(x));

  } // end iteration loop
  cout << endl;

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
			    int LType,
			    TString opt,
			    const GSVDResult* gsvd,
			    const TH2* hA,
			    const TH1* hb,
			    const TH1* hXini)
{
  TMatrixD L = LMatrix(fN, LType);
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
			    TMatrixD& L,
			    TString opt,
			    const GSVDResult* gsvd,
			    const TH2* hA,
			    const TH1* hb,
			    const TH1* hXini)
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
  double x1 = fTrueX1, x2 = fTrueX2;
  if (hA) A = Hist2Matrix(hA);
  if (hb) b = Hist2Vec(hb);
  if (hXini) {
    x1 = hXini->GetXaxis()->GetXmin();
    x2 = hXini->GetXaxis()->GetXmax();
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
  if (W.GetNcols()>0) { // if L has a nontrivial null space
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
  
  for (int k=1; k<=nIterations; k++) {
    printf("PCGLS iteration %d\r", k);
    TVectorD Az = A*z;
    alpha = dq / Az.Norm2Sqr();
    x += alpha*z;
    r -= alpha*Az;
    s  = AT*r;
    LTSolve(q1, L, s);

    // Reorthogonalize q1 to previous q1 vectors
    for (int i=0; i<k; i++) {
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
    if (gsvd) {
      
      if (k==1) {
	fk = alpha*g2;
	Fd = g2 - ElemMult(g2, fk) + beta*g2;
	TMatrixDColumn(result.F, k-1) = fk;
      }
      else {
	fk += alpha*Fd;
	Fd = g2 - ElemMult(g2, fk) + beta*Fd;
	TMatrixDColumn(result.F, k-1) = fk;
      }
      if (k > 2) {
	for (int i=0; i<fN; i++) {
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
      if (gcv < gcvMin) {
	gcvMin = gcv;
	result.kGcv = k;
      }
      
      // Compute A^#, the regularized inverse of A, using these filter
      // factors.
      TMatrixD FCd(fN,fN);
      TMatrixD UT(TMatrixD::kTransposed, gsvd->U);
      for (int i=0; i<fN; i++) {
	double a = gsvd->alpha(i);
	FCd(i,i) = (a > 0.) ? result.F(i,k-1)/a : 0.;
      }
      TMatrixD Ap = gsvd->X * FCd * UT;
      
      // Covariance matrices:
      // Cov(w) = Ap*B*Ap'
      // Cov(x)_ik = xini_i Cov(w)_ik xini_k 
      result.wCov = Ap * TMatrixD(B, TMatrixD::kMultTranspose, Ap);
      for (int i=0; i<fN; i++) {
	for (int j=0; j<fN; j++) {
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
  TH2D* errw = Matrix2Hist(result.WRegErr, Form("errW_CG_%d",id),
			   fTrueX1, fTrueX2,0,nIterations);
  TH2D* errx = Matrix2Hist(result.XRegErr, Form("errX_CG_%d",id),
			   fTrueX1, fTrueX2,0,nIterations);

  for (int j=1; j<=fN; j++) {
    for (int k=1; k<nIterations; k++) {
      result.WRegHist->SetBinError(j,k,errw->GetBinContent(j,k));
      result.XRegHist->SetBinError(j,k,errx->GetBinContent(j,k));
    }
  }

  result.hGcv = 
    result.XRegHist->ProjectionX(Form("cg_%d_iter%d",id,result.kGcv),
				 result.kGcv,result.kGcv);
  
  return result;
}

void 
UnfoldingUtils::LSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y, 
		       const TMatrixD& W, const TMatrixD& T)
{
  // Computes  x = L_p*y
  // where L_p is the A-weighted generalized inverse of L.
  int p = L.GetNrows();
  int n = L.GetNcols();
  int ly = y.GetNoElements();

  if (ly != p)
    gROOT->Warning("LSolve()", 
		   "Input vector length = %d != %d", ly, p);
  if (p==n) {
    TMatrixD Linv(TMatrixD::kInverted, L);
    TVectorD x = Linv * y;
    result = x;
    return;
  }
  
  if (p > n) {
    TMatrixD Ltmp = L;
    TMatrixD Linv = MoorePenroseInverse(Ltmp);
    TVectorD x = Linv * y;
    result = x;
    return;
  }
  else {
    TMatrixD L11(L); L11.ResizeTo(p,p);
    TMatrixD T11(T); T11.ResizeTo(n-p,p);
    TMatrixD L11inv = MoorePenroseInverse(L11);
    TVectorD xhat = L11inv * y;
    TVectorD xhat0 = xhat;
    xhat0.ResizeTo(n); // Append n-p 0's
    TVectorD x = xhat0 - W*T11*xhat; 
    result.ResizeTo(x.GetNoElements());
    result = x;
  }
  return;
 }

void 
UnfoldingUtils::LTSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y)
{
  // Computes  x = (L_p)'*y
  // where L_p is the generalized inverse of L (aka L_A^dagger)
  
  int p = L.GetNrows();
  int n = L.GetNcols();
  int ly = y.GetNoElements();
  
  result.Zero();
  
  if (p == n) {
    TMatrixD LT(L);
    LT.Invert();
    LT.T();
    result = LT * y;
    return;
  }
  if (p > n) {
    TMatrixD LT(L); LT.T();
    result = MoorePenroseInverse(LT) * y;
    return;
  }
  else {

    TVectorD y1(y);
    TMatrixD L11(L);
    y1.ResizeTo(p);
    L11.ResizeTo(p,p);

    // Calculate result
    L11.Invert();
    L11.T();
    result = L11 * y1;

    if (0) {
      Printf("p, n, ly = %d, %d, %d", p, n, ly);
      Printf("y, y1, result: %d, %d, %d", 
	     y.GetNoElements(),
	     y1.GetNoElements(), 
	     result.GetNoElements());
      Printf("L11 (%d x %d)", 
	     L11.GetNrows(), L11.GetNcols());
    }

  }
  return;
 }

TMatrixD 
UnfoldingUtils::MoorePenroseInverse(TMatrixD& A, double tol)
{
  // Compute the Moore-Penrose (pseudo)inverse of A as
  // V*diag(1/sigma_i)*U' (Numerical Recipes eq 2.6.6)
  // Inverse singular values < tol are zeroed.
  int m = A.GetNrows(), n = A.GetNcols();
  int max = m;
  if (n>max) 
    max = n;
  // Get SVD components
  TDecompSVD decomp;
  
  if (m<n) {
    TMatrixD B(A);
    B.ResizeTo(n,n);
    decomp.SetMatrix(B);
  }
  else {
    decomp.SetMatrix(A);
  }
  
  TMatrixD UT = decomp.GetU(); UT.T();
  TVectorD sig = decomp.GetSig();
  TMatrixD V = decomp.GetV();
  TMatrixD Siginv(max,max);
 
  V.ResizeTo(max,max);

  for (int i=0; i<n; i++) {
    if (sig(i) > tol)
      Siginv(i,i) = 1./sig(i);
    else
      Siginv(i,i) = 0.;
  }

  TMatrixD Ainv = (V*Siginv)*UT;
  Ainv.ResizeTo(n,m);

  return Ainv;
}

TGraph*
UnfoldingUtils::LogCurvature(TGraph* g, const TVectorD& tVec, int& kMax)
{
  // Compute signed log curvature of the parametric curve
  // g(t) = (a(t), b(t))
  // where (x(t), y(t)) = (log(a(t)), log(b(t))).
  // Curvature = (x'y'' - y'x'') / (x'^2 + y'^2)^1.5
  // kMax is the index of the maximum point.

  int n = g->GetN();
  TGraph* gk = new TGraph(n-2);

  TVectorD x(n);
  TVectorD y(n);
  for (int i=0; i<n; i++) {
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
  for (int i=0; i<n-1; i++) {
    x1(i) = (x(i+1) - x(i)) / (tVec(i+1) - tVec(i)); 
    y1(i) = (y(i+1) - y(i)) / (tVec(i+1) - tVec(i)); 
  }
  // Second derivatives
  for (int i=0; i<n-2; i++) {
    x2(i) = (x1(i+1) - x1(i)) / (tVec(i+1) - tVec(i)); 
    y2(i) = (y1(i+1) - y1(i)) / (tVec(i+1) - tVec(i)); 
  }

  // Compute curvature vs. t
  double cmax = -1e99;
  for (int i=0; i<n-2; i++) {
    double numer = x1(i)*y2(i) - x2(i)*y1(i);
    double denom = TMath::Power(x1(i)*x1(i) + y1(i)*y1(i), 1.5);
    double c = numer / denom;
    gk->SetPoint(i, tVec(i), c);

    // Find maximum curvature and its array position
    if (c>cmax) {
      cmax = c;
      kMax = i;
    }

  }
  return gk;
}

TVectorD
UnfoldingUtils::Ones(int n)
{
  TVectorD ones(n);
  for (int i=0; i<n; i++)
    ones(i) = 1.0;
  return ones;
}

TVectorD 
UnfoldingUtils::ElemDiv(const TVectorD& x, const TVectorD& y, double div0val)
{
  int nx = x.GetNoElements();
  int ny = y.GetNoElements();
  if (nx != ny) {
    gROOT->Error("ElemDiv()", "mismatch nx=%d, ny=%d", nx, ny);
    gSystem->Exit(-1);
  }
  TVectorD result(nx);
  for (int i=0; i<nx; i++) {
    result(i) = (y(i) > 1e-15) ? x(i) / y(i) : div0val;
  }
  return result;
}

TVectorD 
UnfoldingUtils::ElemMult(const TVectorD& x, const TVectorD& y)
{
  int nx = x.GetNoElements();
  int ny = y.GetNoElements();
  if (nx != ny) {
    gROOT->Error("ElemMult()", "mismatch nx=%d, ny=%d", nx, ny);
    gSystem->Exit(-1);
  }
  TVectorD result(nx);
  for (int i=0; i<nx; i++) {
    result(i) = x(i)*y(i);
  }
  return result;
}

TMatrixD 
UnfoldingUtils::DivColsByVector(const TMatrixD& M, const TVectorD& v,
				bool makeZeroIfNaN)
{
  // Divide M columns by v elementwise: R(i,j) = M(i,j) / v(i).
  // I.e. each row i is scaled by 1/v(i). 
  TMatrixD R(M);
  int m = R.GetNrows(), n = R.GetNcols();

  if (v.GetNoElements() != m)
    Error("UnfoldingUtils::DivColsByVector()", 
	  "nrows %d != vector size %d", m, v.GetNoElements());
  
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      if (v(i) != 0) 
	R(i,j) /= v(i);
      else if (makeZeroIfNaN) 
	R(i,j) = 0;
    }
  }
  return R;
}

TMatrixD 
UnfoldingUtils::MultRowsByVector(const TMatrixD& M, const TVectorD& v)
{
  // Multiply M rows by v elementwise: R(i,j) = M(i,j) * v(j).
  // I.e. each column j is scaled by v(j). 
  TMatrixD R(M);
  int m = R.GetNrows(), n = R.GetNcols();

  if (v.GetNoElements() != n)
    Error("UnfoldingUtils::MultRowsByVector()", 
	  "ncols %d != vector size %d", n, v.GetNoElements());
  
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
	R(i,j) *= v(j);
    }
  }
  return R;
}


TH1D* 
UnfoldingUtils::XHist(TVectorD& x, TString base, int k, 
		      double xMin, double xMax,
		      double normto, TString /*opt*/)
{
  // base should have a form like "CGLS"+"id"
  // k is the iteration step
  // The factor normto could be hb->Integral();
  const char* name  = Form("h%s_%d", base.Data(), k);
  const char* title = Form("%s (%d iterations)", base.Data(), k);
  TH1D* hx = Vec2Hist(x, xMin, xMax, name, title);

  if (normto != 0.0) {
    double norm = hx->Integral() ? normto/hx->Integral() : 0;
    hx->Scale(norm);
  }
  // TODO assign proper stat uncertainty
  return hx;
}

TMatrixD 
UnfoldingUtils::Toeplitz(int m1, int n1, double col[], double row[])
{
  // Return a Toeplitz matrix T constructed from col[] (size m+1) and
  // row[] (size n+1). T has dimensions (m+1) x (n+1). T is assigned
  // as
  //
  //                   c[0] = r[0], i==j (diagonal)
  //    t[i][j] =      c[i-j],      i>j  (lower left)
  //                   r[j-i],      i<j  (upper right)
  //
  // To create arrays to pass in, follow this pattern:
  // c[m+1] = {t[0], t[-1], t[-2], ..., t[-m]} (i.e. t[-i] = c[i])
  // r[n+1] = {t[0], t[1],  t[2],  ..., t[n]}  (i.e. t[j]  = r[j])

  if (col[0] != row[0]) {
    Warning("UnfoldingUtils::Toeplitz()", 
	    "col[0] (%f) != row[0] (%f). Using col[0].", 
	    col[0], row[0]);
  }
  
  TMatrixD T(m1, n1);
  T(0, 0) = col[0];
  for (int i=0; i<m1; i++) {
    for (int j=0; j<n1; j++) {
      if     (i==j)  T(i, j) = col[0];
      else if (i>j)  T(i, j) = col[i-j];
      else if (i<j)  T(i, j) = row[j-i];
    }
  }
  return T;
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

  if (kind == kUnitMatrix) { // n x n identity matrix. p = n. W=0.
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    return Toeplitz(n,n,col,row);
  }
  if (kind == k1DerivNoBC) { // 1st deriv, no BC assumptions. p = n-1. W=const.
    double col[n-1];
    for (int i=0; i<n-1; i++) col[i]=0.0;
    row[0] = col[0] = -1;
    row[1] = 1;
    return Toeplitz(n-1,n,col,row);
  }
  if (kind == k2DerivNoBC) { // 2nd deriv, no BC assumptions. p = n-2. W=const, linear.
    double col[n-2];
    for (int i=0; i<n-2; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    row[1] = -2;
    row[2] = 1;
    return Toeplitz(n-2,n,col,row);
  }
  if (kind == k1DerivBC0) { // 1st deriv, BC=0 on L and R. p = n+1. W=0.
    double col[n+1];
    for (int i=0; i<n+1; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    col[1] = -1;
    return Toeplitz(n+1,n,col,row);
  }
  if (kind == k2DerivBC0) { // 2nd deriv, BC=0 on L and R. p = n. W=0.
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = -2+eps;
    row[1] = col[1] = 1;
    return Toeplitz(n,n,col,row);
  }
  if (kind == k1DerivBCR) { // 1st deriv, reflect at L,R. p = n-1. W=const. Same as k1DerivNoBC
    double col[n-1];
    for (int i=0; i<n-1; i++) col[i]=0.0;
    row[0] = col[0] = -1;
    row[1] = 1;
    return Toeplitz(n-1,n,col,row);
  }
  if (kind == k2DerivBCR) { // 2nd deriv, reflect at L,R. p = n. W=const.
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

int
UnfoldingUtils::Rank(TMatrixD& A)
{
  // Not very efficient, could replace with a rank-revealing QR
  // algorithm.
  TMatrixD N = Null(A);
  int nullity = N.GetNcols();
  int rank = TMath::Min(A.GetNrows(), A.GetNcols()) - nullity;
  return rank;
}

TMatrixD 
UnfoldingUtils::Null(TMatrixD& A)
{
  // Return matrix whose columns form the orthonormal basis for the
  // nullspace of A. Perform SVD on A so that A = USV'. The columns of
  // V whose corresponding singular values are "zero" are the null basis
  // vectors.
  // TODO: pass in tolerance instead of hardcoding 1e-16
  int m = A.GetNrows(), n = A.GetNcols();

  // Result
  TMatrixD W(0, 0);

  // Get SVD components
  TDecompSVD decomp;

  if (m<n) {
    TMatrixD B = A;
    B.ResizeTo(n,n);
    decomp.SetMatrix(B);
  }
  else {
    decomp.SetMatrix(A);
  }

  TVectorD sig = decomp.GetSig();
  TMatrixD V = decomp.GetV();
  
  if (0) {  
    sig.Print();
    V.Print();
  }
  
  int ndim = 0;
  for (int i=0; i<sig.GetNoElements(); i++) {
    if (sig(i) < 1e-16) {
      ndim++;
      W.ResizeTo(n, ndim);
      TMatrixDColumn(W,ndim-1) = TMatrixDColumn(V, i);
    }
  }
  
    return W;
}

void 
UnfoldingUtils::SetTH1Props(TH1* h,
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

TH2D* 
UnfoldingUtils::BandedDiagonalMatrix(TH1* hDpt, 
				     const int nMeas, 
				     const int nTrue, 
				     double xm1, 
				     double xm2, 
				     double xt1, 
				     double xt2)
{
  // Make a bkg matrix by repeating dpt bins along diagonals.
  // zerobin holds the (0,0) element, = firstCol[0] = firstRow[0]. 
  // nMeas must be <= the number of positive pt bins.
  // nTrue must be <= the number of negative pt bins.
   
  int zerobin = hDpt->FindBin(1e-6); // index of 1st positive bin
  double dptIntegral = hDpt->Integral();
  double wdpt = hDpt->GetBinWidth(1); 
  double firstRow[9999] = {0};
  double firstCol[9999] = {0};
  int ip=0;
  
  // Check for acceptable setup
  int nNegBins = zerobin-1;          // binning along true TH2 axis
  int nPosBins = hDpt->GetNbinsX() - nNegBins;  // "  meas TH2 axis
  if (nNegBins < 1)
    Error("UnfoldingUtils::BandedDiagonalMatrix()", 
	  "nNegBins < 1 (%d)", nNegBins);
  if (nMeas > nPosBins)
    Error("UnfoldingUtils::BandedDiagonalMatrix()", 
	  "nMeas (%d) > nPosBins (%d)", nMeas, nPosBins);
  if (nTrue >= zerobin)
    Error("UnfoldingUtils::BandedDiagonalMatrix()", 
	  "nTrue (%d) > nNegBins (%d)", nTrue, nNegBins);
  
  ip=0;
  for (int k=zerobin; k<=zerobin+nMeas; k++) {
    // Set content of leftmost column (positive pT bins, TMatrixD true axis).
    firstCol[ip++] = hDpt->GetBinContent(k) / dptIntegral;
  }
  ip=0;
  for (int k=zerobin; k<=zerobin+nTrue; k++) {
    // Set content of first row (negative pT bins, TMatrixD meas axis).
    firstRow[ip] = hDpt->GetBinContent(zerobin-ip) / dptIntegral;
    ip++;
  }

  TMatrixD Abkg = Toeplitz(nMeas,nTrue,firstCol,firstRow);
  double x1=0, y1=0, x2=wdpt*nMeas, y2=wdpt*nTrue;
  if (xm2 != xm1) {
    x1 = xm1;
    x2 = xm2;
  }
  if (xt2 != xt1) {
    y1 = xt1;
    y2 = xt2;
  }
  
  TH2D* h = Matrix2Hist(Abkg, "banded_diag",x1,x2,y1,y2);
  h->SetTitle("A_{bkg} from #deltap_{T};"
	      "Reconstructed p_{T} (GeV/c);"
	      "Generated p_{T} (GeV/c)");
  return h;
}

TH2*
UnfoldingUtils::TH2Product(TH2* hA, TH2* hB, TString name)
{
  // Return the matrix product hA*hB as a new TH2
  int nA = hA->GetNbinsY(), mB = hB->GetNbinsX();

  if (nA != mB) {
    Error("UnfoldingUtils::TH2Product()", 
	  "Number of A cols (%d) != B rows (%d)", nA, mB);
    return 0;
  }
  
  double xa1 = hA->GetXaxis()->GetXmin(), xa2 = hA->GetXaxis()->GetXmax();
  double yb1 = hB->GetYaxis()->GetXmin(), yb2 = hB->GetYaxis()->GetXmax();

  TMatrixD A = Hist2Matrix(hA);
  TMatrixD B = Hist2Matrix(hB);
  TMatrixD AB(A*B);
  return Matrix2Hist(AB, name, xa1, xa2, yb1, yb2);
}

TMatrixD
UnfoldingUtils::OuterProduct(TVectorD a, TVectorD b)
{
  // Return a*b'
  int na = a.GetNrows();
  int nb = b.GetNrows();
  TMatrixD M(na,nb);
  for (int i=0; i<na; i++) {
    for (int j=0; j<nb; j++) {
      M(i,j) = a(i)*b(j);
    }
  }
  return M;
}

TH2D*
UnfoldingUtils::TH2Sub(TH2* h, int bx1, int bx2, int by1, int by2, TString name)
{
  // Return a sub-range of h (bx1..bx2) x (by1..by2), bx2,by2 included.
  double xw = h->GetXaxis()->GetBinWidth(1);
  double yw = h->GetYaxis()->GetBinWidth(1);
  double x1 = (bx1-1)*xw, x2 = bx2*xw;
  double y1 = (by1-1)*yw, y2 = by2*yw;
  int nx = bx2-bx1+1, ny = by2-by1+1;

  if (bx1<1 || bx2>h->GetNbinsX()) {
    Error("UnfoldingUtils::TH2Sub()", 
	  "Requested x bins %d,%d out of range (1,%d)", bx1,bx2,h->GetNbinsX());
    return 0;
  }
  if (by1<1 || by2>h->GetNbinsY()) {
    Error("UnfoldingUtils::TH2Sub()", 
	  "Requested y bins %d,%d out of range (1,%d)", by1,by2,h->GetNbinsY());
    return 0;
  }

  TH2D* hSub = new TH2D(name.Data(),name.Data(),nx,x1,x2,ny,y1,y2);

  int ipx=1, ipy=1;
  for (int ix=bx1; ix<=bx2; ix++) {
    for (int iy=by1; iy<=by2; iy++) {
      double val = h->GetBinContent(ix, iy);
      hSub->SetBinContent(ipx, ipy, val);
      ipy++;
    }
    ipy=1;
    ipx++;
  }
  
  return hSub;
}

TMatrixD
UnfoldingUtils::DerivativeMatrix(int n, int d)
{
  // Return a d-th derivative operator matrix.
  // n is the number of columns (as usual).
  // Default size is n-d x n.
  // Unit matrix is returned if d=0.

  int nd = n-d;
  TMatrixD L(nd, n);
  TVectorD c(d+1); // Use to populate L
  c.Zero();

  // Assign c  
  if (d==0)
    c(0) = 1;
  else {
    c(0) = -1;
    c(1) = 1;
  }
  for (int i=2; i<=d; i++) {
    TVectorD a(d+2), b(c);  
    a.SetSub(1, c);
    a.ResizeTo(d+1);
    b(d) = 0;
    c = a-b;
  }

  for (int i=0; i<nd; i++) {
    for (int j=0; j<d+1; j++) {
      L(i,j+i) = c(j);
    }
  }
  
  return L;
}

void 
UnfoldingUtils::NormalizeXSum(TH2* hA, TH1* hN)
{
  // Normalize x-rows of hA so each row sums to 1.0 (default), or
  // optionally to the value of the jth bin of hN.
  
  int nx = hA->GetNbinsX();
  int ny = hA->GetNbinsY();
  
  if (hN)
    if (ny != hN->GetNbinsX())
      Error("NormalizeXSum()", 
	    "ny=%d != %d in hN", nx, hN->GetNbinsX());
  
  // xsum(j) contains sum of x cells in row j
  TVectorD xsum(ny);
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      xsum(j) += hA->GetBinContent(i+1,j+1);
    }
  }
  
  // Change bin contents of hA to normalized value, which is 1.0 if hN
  // is not passed in.
  for (int j=0; j<ny; j++) {
    double a = hN ? hN->GetBinContent(j+1) : 1.; 
    double w = (xsum(j) != 0.) ? a/xsum(j) : a;
    for (int i=0; i<nx; i++) {
      double val = w*hA->GetBinContent(i+1,j+1);
      hA->SetBinContent(i+1,j+1, val);
    }
  }
}

TGraph*
UnfoldingUtils::ResidualNorm(TObjArray* hists, double stepSize)
{
  static int id=0; id++;

  // Check for Ax-b histos in hists
  TObjArray* subList = new TObjArray();
  for (int i=0; i<hists->GetEntries(); i++) {
    TObject* obj = hists->At(i);
    TString name = obj->GetName();
    TString cl   = obj->ClassName();
    if (name.Contains("hr") && cl.Contains("TH1")) {
      subList->Add(obj);
      if (0) 
	Info("UnfoldingUtils::ResidualNorm()", "Added %s",name.Data());
    }
  }
  int np = subList->GetEntries();
  TGraph* g = new TGraph(np);
  for (int i=0; i<np; i++) {
    TH1* h = (TH1*)subList->At(i);
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

TH2D*
UnfoldingUtils::UnfoldCovMatrix(int nTrials, 
				int algo, 
				double regPar,
				TString opt,
				TObjArray* bHists)
{
  // Propagate the measured data cov matrix to the unfolded result.
  static int id=0; id++;
  
  if (bHists && nTrials > bHists->GetEntries()) {
    Error("UnfoldingUtils::UnfoldCovMatrix()",
	  "nTrials (%d) > number of trial data sets (%d)",
	  nTrials, bHists->GetEntries());
    return 0;
  }
  
  // seed = 0 signals to use the TUUID identifier
  int seed = 0;
  TRandom3 rand(seed);
  
  // Unfolded result
  TH1D* hUnf = 0;
  
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
  TH2D* hcov = new TH2D(Form("hcov%d",id),title.Data(),
			fN,fTrueX1,fTrueX2,fN,fTrueX1,fTrueX2); 

  // Gaussian perturbation
  TVectorD g(fM);

  // Measured data + MC perturbation for one pseudo-experiment
  TH1D* hbTrial = (TH1D*)fHistMeas->Clone(Form("hbTrial%d",id));

  // Mean result of unfolding MC trials E[x]
  TH1D* hxMean = new TH1D(Form("hxMean%d",id),"mean",fN,fTrueX1,fTrueX2); 
  
  // Loop 1: Compute E[x] as hxMean
  for (int t=0; t<nTrials; t++) { 
 
    cout << Form("UnfoldCovMatrix(): computing mean result %d\r",t) << flush;

    if (bHists) {
      fHistMeas = (TH1D*)bHists->At(t);
      fVecb = Hist2Vec(fHistMeas);
      fVecbErr = Hist2Vec(fHistMeas, "unc");
      fMatATilde = DivColsByVector(fMatA, fVecbErr);
      fVecbTilde = ElemDiv(fVecb, fVecbErr);
      
      /////
    }
    else {
      for (int i=0; i<fM; i++) {
	g(i) = rand.Gaus(0,1);
      }
      
      g *= B;//UT;
      for (int i=0; i<fM; i++) {
	hbTrial->SetBinContent(i+1, fVecb(i) + g(i));
	hbTrial->SetBinError(i+1, fHistMeas->GetBinError(i+1));
	fVecb(i) = hbTrial->GetBinContent(i+1);
	fVecbErr(i) = hbTrial->GetBinError(i+1);
      }
    }

    switch (algo) {
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
      // 	h = UnfoldChi2Min(); 
      // 	break;
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
  for (int t=0; t<nTrials; t++) {  

    cout << Form("UnfoldCovMatrix(): computing covariance %d\r",t) << flush;

    if (bHists) {
      fHistMeas = (TH1D*)bHists->At(t);
      fVecb = Hist2Vec(fHistMeas);
      fVecbErr = Hist2Vec(fHistMeas, "unc");
      fMatATilde = DivColsByVector(fMatA, fVecbErr);
      fVecbTilde = ElemDiv(fVecb, fVecbErr);
    }
    else {
      for (int i=0; i<fM; i++) {
	g(i) = rand.Gaus(0,1);
      }
      g *= B;//UT;
      for (int i=0; i<fM; i++) {
	hbTrial->SetBinContent(i+1, fVecb(i) + g(i));
	hbTrial->SetBinError(i+1, fHistMeas->GetBinError(i+1));
      }
    }
    switch (algo) {
    case kGSVDAlgo: 
      hUnf = GSVDAnalysis(fMatL, regPar, 0, 0, opt)->xregHist;
      break;
    case kSVDAlgo: 
      hUnf = UnfoldSVD(regPar, 0, opt, 0, 0, 0); 
      break;
    }
    
    double dxj=0, dxk=0, cjk=0;
    for (int j=1; j<=fN; j++) {
      for (int k=1; k<=fN; k++) {
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

QRDecompResult 
UnfoldingUtils::QRDecomp(TMatrixD& A)
{
  // Compute QR decomposition of A using Householder transformations
  QRDecompResult qr;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD R(A);
  TMatrixD Qj(m,m);

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++) {
    TVectorD col = TMatrixDColumn(R,j);
    TVectorD x = col.GetSub(j,m-1);
 
    int sign = (col(j)<0.)? -1. : 1.;
    double alpha = sign*TMath::Sqrt(x*x);
    TVectorD u(x);
    u(0) += alpha;

    // Compute Householder vector v and matrix H
    double unorm = TMath::Sqrt(u*u);
    TVectorD v(u); v *= (unorm==0)? 0. : 1./unorm;  
    TMatrixD H = OuterProduct(v,v);
    H *= 2;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;
    
    // Full-dimension Householder matrix
    Qj.UnitMatrix();
    Qj.SetSub(j,j,H);

    // Update Q and R    
    Q = Q*Qj;    
    R = Qj*R;
  }

  bool requirePositivePivotEntries = true;
  if (requirePositivePivotEntries) {
    int r = TMath::Min(m,n);
    for (int i=0; i<r; i++) {
      if(R(i,i)<0) {
	TMatrixDRow(R,i) *= -1;
	TMatrixDColumn(Q,i) *= -1;
      }
    }
  }

  // Store results
  qr.Q.ResizeTo(Q);
  qr.R.ResizeTo(R);
  qr.Q = Q;
  qr.R = R;
  
  return qr;
}

QRDecompResult 
UnfoldingUtils::QLDecomp(TMatrixD& A)
{
  // Compute QL decomposition of A using Householder transformations
  QRDecompResult ql;
  int m = A.GetNrows();
  int n = A.GetNcols();
  TMatrixD Q(m,m); Q.UnitMatrix();
  TMatrixD L(A);
  TMatrixD Qj(m,m);

  // For sign manipulation
  TMatrixD U1(m,m);
  U1.UnitMatrix();
  TMatrixD L1(n,n);
  L1.UnitMatrix();

  int nIter = TMath::Min(m-1, n);
  for (int j=0; j<nIter; j++) {
    TVectorD col = TMatrixDColumn(L,n-j-1);
    TVectorD x = col.GetSub(0,nIter-j);

    int sign = (x(nIter-j)<0.)? -1. : 1.;
    double alpha = sign*TMath::Sqrt(x*x);
    TVectorD u(x);
    u(nIter-j) += alpha;

    // Compute Householder vector v and matrix H
    double unorm = TMath::Sqrt(u*u);
    TVectorD v(u); v *= (unorm==0)? 0. : 1./unorm;  
    TMatrixD H = OuterProduct(v,v);
    H *= 2;
    TMatrixD I(H); I.UnitMatrix();
    H = I - H;

    // Full-dimension Householder matrix
    Qj.UnitMatrix();
    Qj.SetSub(0,0,H);
    
    // Update Q and L
    Q = Q*Qj;    
    L = Qj*L;
  }

  // cout<<"Q before"; Q.Print();
  // cout<<"L before"; L.Print();
  
  bool requirePositivePivotEntries = true;
  if (requirePositivePivotEntries) {
    int r = TMath::Min(m,n);
    int d = m-n;
    for (int i=0; i<r; i++) {
      int row=i,col=i;
      if (m<n)
	col -= d;
      else if (m>n)
	row += 1;
      
      //      Printf("i= %d, L(%d,%d) = %g", i, row, col, L(row,col));
      if(L(row,col)<0) {
	TMatrixDRow(L,row) *= -1;
	TMatrixDColumn(Q,row) *= -1;
      }
    }
  }
  
  // Store results
  ql.Q.ResizeTo(Q);
  ql.R.ResizeTo(L);
  ql.Q = Q;
  ql.R = L;
  
  return ql;
}


CSDecompResult 
UnfoldingUtils::CSDecomp(TMatrixD& Q1, TMatrixD& Q2)
{
  // ***************************************
  // Q1 shorter or equal to Q2 (m <= p case)
  // ***************************************
  bool debug = false;
  int m,p,l,q1,q2,r,n;
  m  = Q1.GetNrows();
  p  = Q2.GetNrows();
  l  = Q1.GetNcols();
  r  = 0;

  if (m > p) {
    Error("UnfoldingUtils::CSDecomp()",
	  "Q1 rows (%d) <= Q2 (%d) rows required.\nExiting.",m,p);
        gSystem->Exit(-1);
  }

  TMatrixD C(m,l);
  TMatrixD S(p,l);
  CSDecompResult csd;
  TVectorD alpha(l);
  TVectorD beta(l);

  // 1.
  q1 = TMath::Min(m,l); // 3
  q2 = TMath::Min(p,l); // 4

  // 2. SVD of Q2: Q2 = VSZ'
  TDecompSVD svdQ2(Q2);
  TMatrixD V    = svdQ2.GetU();    // p x p
  TMatrixD Z    = svdQ2.GetV();    // l x l
  beta = svdQ2.GetSig();           // l

  // 3-5. Re-order V, S, Z cols
  ReverseColumns(V,0,q2-1);
  ReverseColumns(Z);
  ReverseVector(beta);

  // 6.
  for (int i=0; i<l-q2; i++) {
    alpha(i) = 1.0;
    beta(i)  = 0.0; 
  }

  if (debug) {  
    cout << "V: ";  V.Print();
    //    cout << "S: ";  S.Print();
    cout << "Z: ";  Z.Print();
    cout << "beta:"; beta.Print(); // non-decreasing
  }
  
  // 7.
  // Find r where beta(r) <= 1/sqrt(2) < beta(r+1).  C++ problem:
  // index != dimension! Need two variables, r and rdim, to resolve
  // the ambiguity when r=0. rdim is the number of betas below 0.707,
  // r is the index.
  double thr = 1./TMath::Sqrt(2.);
  int rdim = 0;
  r = -1;
  for (int i=0; i<m-1; i++) {
    if (beta(i) <= thr && beta(i+1) > thr) {
      r = i;
      break;
    }
  }
  if (r == -1) {
    r = 0;
  }
  else rdim = r+1;

  if (debug) {  
    Printf("r = %d, rdim = %d.",r,rdim);
  }
  
  // 8.
  TMatrixD T = Q1*Z; // (m x l)

  // 9.
  // QR decomp of T: T = UR
  QRDecompResult qrT = QRDecomp(T);
  TMatrixD U = qrT.Q;
  TMatrixD R = qrT.R;
  
  if (debug) {  
    cout << "T: ";  T.Print();
    cout << "U and R: ";
    U.Print();
    R.Print();
  }

  // Get R2 and R3 from R
  TMatrixD R2 = R.GetSub(l-q2,r,l-q2,r);
  TMatrixD R3 = R.GetSub(rdim,q1-1,rdim,l-1);
  int r3r = R3.GetNrows();
  int r3c = R3.GetNcols();
  if (r3r < r3c)
    R3.ResizeTo(r3c, r3c);

  if (debug) {  
    cout << "R2: ";  R2.Print();
    cout << "R3: ";  R3.Print();
  }
  
  // 10.
  // Compute SVD of R3: R3 = Ur*Cr*Zr'
  TDecompSVD svd2(R3);
  TMatrixD Ur = svd2.GetU();
  TMatrixD Zr = svd2.GetV();
  TVectorD a3 = svd2.GetSig();

  for (int i=0; i<a3.GetNrows(); i++)
    alpha(rdim+i) = a3(i);

  // 11.
  for (int i=q1; i<l; i++) {
    alpha(i) = 0.0;
    beta(i)  = 1.0;
  }

  // 12.
  for (int i=l-q2; i<rdim; i++) {
    alpha(i) = R2(i,i);
  }

  if (debug) {  
    cout << "alpha: ";  alpha.Print();
  }
  
  // 13.
  // Form final U matrix
  // First resize U to undo TDecompSVD-required modification
  if (r3r < r3c) {
    Ur.ResizeTo(r3r,r3r);
  }

  TMatrixD Ur1(U);
  Ur1.UnitMatrix();
  Ur1.SetSub(rdim,rdim,Ur);

  if (debug) {
    cout << "Ur: ";  Ur.Print();
    cout << "Ur1: ";  Ur1.Print();
  }

  // 14.
  // Form final Z matrix
  TMatrixD Zrt(Zr);
  Zrt.ResizeTo(q2-rdim,q2-rdim);

  TMatrixD Zr1(Z);
  Zr1.UnitMatrix();  
  Zr1.SetSub(rdim,rdim,Zr);

  if (debug) {
    cout << "Zr1: ";  Zr1.Print();
  }
  
  U = U*Ur1;
  Z = Z*Zr1;

  // 15.
  TMatrixD St(Zrt);
  St.Zero();
  for (int i=0; i<St.GetNcols(); i++)
    St(i,i) = beta(i+rdim);

  TMatrixD W = St*Zrt;
  
  // 16.
  QRDecompResult qrw = QRDecomp(W);
  TMatrixD Qw = qrw.Q;
  n = TMath::Min(rdim,l-q2);

  // 17.
  TMatrixD Vpost(V);
  Vpost.UnitMatrix();
  Vpost.SetSub(rdim-n, rdim-n, Qw);
  V = V*Vpost;
  
  // Construct C from alpha
  for (int i=0; i<q1; i++)
    C(i,i) = alpha(i);

  // And S from beta
  for (int i=0; i<q2; i++)
    S(i,i) = beta(i);

  if (debug) {
    cout << "St: ";  St.Print();
    cout << "Zrt: ";  Zrt.Print();
    cout << "W: ";  W.Print();
    cout << "Qw: ";  Qw.Print();
    cout << "Vpost: ";  Vpost.Print();

    cout << "C: ";  C.Print();
    cout << "S: ";  S.Print();
    TMatrixD CC(C,TMatrixD::kTransposeMult,C);
    TMatrixD SS(S,TMatrixD::kTransposeMult,S);
    TMatrixD one = CC + SS;
    cout << "C'C + S'S: ";  one.Print();

    TMatrixD G(S,TMatrixD::kMultTranspose,Z);
    TMatrixD Ginv = MoorePenroseInverse(G);
    TMatrixD Vtest = Q2*Ginv;
    cout << "Vtest: ";  Vtest.Print();
    cout << "V: ";  V.Print();

    TMatrixD zz = Q2 - Vtest*G;
    cout << "zz: ";  zz.Print();
  }
  
  csd.C.ResizeTo(C);
  csd.S.ResizeTo(S);
  csd.U.ResizeTo(U);
  csd.V.ResizeTo(V);
  csd.Z.ResizeTo(Z);
  csd.alpha.ResizeTo(l);
  csd.beta.ResizeTo(l);

  csd.C = C;
  csd.S = S;
  csd.U = U;
  csd.V = V;
  csd.Z = Z;
  csd.alpha = alpha;
  csd.beta = beta;

  if (debug)
    Printf("m=%d, p=%d, l=%d, q1=%d, q2=%d, r=%d, n=%d",  m,p,l,q1,q2,r,n);
  return csd;
}

CSDecompResult 
UnfoldingUtils::CSDecompQ1Taller(TMatrixD& Q1, TMatrixD& Q2)
{
  // ***************************************
  // m > p case
  // ***************************************
  bool debug = false;
  int m,p,l,q1,q2,r,rdim;
  m  = Q1.GetNrows();
  p  = Q2.GetNrows();
  l  = Q1.GetNcols();
  rdim = 0;  // number of alpha values < 0.707
  r  = -1;   // index of 1st alpha < 0.707
  TMatrixD C(m,l);
  TMatrixD S(p,l);
  CSDecompResult cs;
  TVectorD alpha(l);
  TVectorD beta(l);

  // 1.
  q1 = TMath::Min(m,l);
  q2 = TMath::Min(p,l);

  // 2.
  // SVD of Q1: UCZ'
  if (fVerbosity)
    Info("UnfoldingUtils::CSDecompQ1Taller()",
	 "Computing SVD of Q1 (%d x %d)",m,l);
  TDecompSVD svdQ1(Q1);
  TMatrixD U     = svdQ1.GetU();    // m x m
  TMatrixD Z     = svdQ1.GetV();    // l x l
  alpha          = svdQ1.GetSig();  // l

  // 3.
  for (int i = 0; i<l; i++) {
    beta(i)  = 1.; 
    if (i>=q1) alpha(i) = 0.;
  }
  
  if (fVerbosity)
    Info("UnfoldingUtils::CSDecompQ1Taller()",
	 "m=%d, p=%d, l=%d, q1=%d, q2=%d",  m,p,l,q1,q2);
  if (debug) {
    Printf("\nSVD: Q1 = U*diag(alpha)*Z\'");
    cout << "U: ";  U.Print();
    cout << "alpha: ";  alpha.Print();
    cout << "Z: ";  Z.Print();
  }
  
  // 4.
  // Find r where alpha(r) >= 1/sqrt(2) > alpha(r+1)
  // r is the index, rdim is the # of values. 
  double thr = 1./TMath::Sqrt(2.);
  for (int i=0; i<l-1; i++) {
    // if (alpha(i) >= thr && alpha(i+1) < thr) {
    //   r = i;
    if (alpha(i) >= thr)
      r = i;
    if (alpha(i+1) < thr)
      break;
  }
  if (r == -1) {
    r = 0;
  }
  else rdim = r+1;
  
  if (debug) 
    Printf ("r = %d, rdim = %d",r,rdim);

  // 5.
  TMatrixD T = Q2*Z;

  // 6.
  // QL decomp of T: T = VL
  QRDecompResult vl = QLDecomp(T);

  TMatrixD V = vl.Q;
  TMatrixD L = vl.R;

  if (debug) {
    cout << "T = Q2*Z = V*L: ";  T.Print();
    cout << "V: ";  V.Print();
    cout << "L (before permutation): ";  L.Print();
    TMatrixD Tcheck(T);
    Tcheck -= V*L;
    Printf("T - V*L sum: %g",Tcheck.Sum());
  }
  
  // Create permutation matrix Pi; L = Pi*L.
  TMatrixD Pi(p,p);
  TMatrixD Iq2(q2,q2); Iq2.UnitMatrix();
  Pi.SetSub(0,p-q2,Iq2);
  if (p>q2) {
    int d = p-q2;
    TMatrixD Id(d,d); Id.UnitMatrix();
    Pi.SetSub(p-1,0,Id);
  }
  L = Pi*L;

  // And its inverse (for later)
  TMatrixD PiInv(TMatrixD::kInverted, Pi);

  if (debug) {
    cout << "Pi: ";  Pi.Print();
    cout << "PiInv: ";  PiInv.Print();
    cout << "Pi*L: ";  L.Print();
  }

  // Get [ L_11 L_12 ], call it L1
  //  int lastrow = (r>0)? r-1 : 0;
  int l1r = (r>0)? p-q2+r-1 : p-q2;
  int l1c = (r>0)? l-q2+r-1 : l-q2;

  if (l1r > L.GetNrows()-1)
    l1r = L.GetNrows()-1;
  if (l1c > L.GetNcols()-1)
    l1c = L.GetNcols()-1;
  
  if (debug)
    Printf("Getting L1 (%d x %d) from L (%d x %d)", l1r,l1c,L.GetNrows(),L.GetNcols());

  TMatrixD L1 = L.GetSub(p-q2, l1r, p-q2, l1c);

  // Get L_23 := L2
  int row1 = p-q2+r; 
  if (row1 > L.GetNrows()-1) 
    row1 = L.GetNrows()-1;
  if (row1 < 0) row1 = 0;
  int col1 = l-q2+r;
  if (col1 > L.GetNcols()-1) 
    col1 = L.GetNcols()-1;
  if (col1 < 0) col1 = 0;
  
  if (L1.GetNrows() < L1.GetNcols()) { // So TDecompSVD works
    if (debug)
      Printf("Resizing L1 (%d x %d) --> (%d x %d)",
	     L1.GetNrows(), L1.GetNcols(), L1.GetNcols(), L1.GetNcols());
    L1.ResizeTo(L1.GetNcols(), L1.GetNcols());
  }
  
  if (debug) {
    cout << "L1 = [L_11 L_12] = Vl*Sl*Zl\': ";  L1.Print();
    Printf("Getting L2... row %d to %d, col %d to %d",
	   row1, L.GetNrows()-1, col1, L.GetNcols()-1);
  }
  
  TMatrixD L2 = L.GetSub(row1, L.GetNrows()-1, col1, L.GetNcols()-1);
  
  if (debug) {
    cout << "L2: ";  L2.Print();
  }
  
  // 7.
  if (fVerbosity)
    Info("UnfoldingUtils::CSDecompQ1Taller()",
	 "Computing SVD of L1 (%d x %d)",L1.GetNrows(),L1.GetNcols());
  TDecompSVD svdl(L1);
  TMatrixD Vlbig = svdl.GetU();
  
  int rmax = TMath::Min(r-1, V.GetNrows()-1);

  if (debug)
    Printf("Getting Vl from Vlbig (%d x %d)... row 0 to %d, col 0 to %d",
	   Vlbig.GetNrows(), Vlbig.GetNcols(), rmax,rmax);

  TMatrixD Vl = Vlbig.GetSub(0,rmax,0,rmax);
  TVectorD bl = svdl.GetSig();
  TMatrixD Zl = svdl.GetV(); 

  // 8-10.
  ReverseVector(bl);
  ReverseColumns(Vl);
  ReverseColumns(Zl);

  // 11.
  int imax = TMath::Min(l, bl.GetNrows());
  for (int i=0; i<imax; i++) 
    beta(i) = bl(i);
  int nl2 = L2.GetNrows();
  for (int i=0; i<nl2; i++) {
    if (bl.GetNrows()+i < l)
      beta(bl.GetNrows()+i) = L2(i,i);
  }
  for (int i=0; i<l-q2; i++) {
    alpha(i) = 1.;
    beta(i)  = 0.; 
  }

  // 12. Create V~ to post-multiply V
  TMatrixD Vpost(V);
  Vpost.UnitMatrix();
  Vpost.SetSub(p-q2,p-q2,Vl);

  // 13.
  if (debug)
    Printf("V = V (%d x %d) * Vpost( %d x %d) * PiInv(%d x %d)",
	   V.GetNrows(), V.GetNcols(), 
	   Vpost.GetNrows(), Vpost.GetNcols(), 
	   PiInv.GetNrows(), PiInv.GetNcols());
  V = V*Vpost*PiInv;
  
  if (debug) {
    cout << "bl (after reversing):";  bl.Print();
    cout << "beta:";  beta.Print();
    cout << "Vl: ";  Vl.Print();
    cout << "Vpost: ";  Vpost.Print();
    cout << "V (= V*Vpost, final): ";  V.Print();
    cout << "Zl: ";  Zl.Print();
    cout << "Z: ";  Z.Print();
  }
  
  // 14.
  TMatrixD Zpost(Z);
  Zpost.UnitMatrix();
  Zpost.SetSub(0,0,Zl);
  Z = Z*Zpost;

  if (debug) {
    cout << "Z (final): "; Z.Print();
  }

  // 15. W = S~ * Zl
  TMatrixD W(Zl);
  W.Zero();
  for (int i=0; i<Zl.GetNcols(); i++) 
    W(i,i) = alpha(i);
  
  if (debug) {
    cout << "S~: ";  W.Print();
  }

  W = W*Zl;
  if (debug) {
    cout << "W: ";  W.Print();
  }
  
  if (fVerbosity)
    Info("UnfoldingUtils::CSDecompQ1Taller()",
	 "Computing QRDecomp(W) (%d x %d)",W.GetNrows(),W.GetNcols());
  QRDecompResult qrw = QRDecomp(W);
  TMatrixD Qw = qrw.Q;
  int ndiff =  U.GetNcols() - Qw.GetNrows(); 
  if (ndiff > 0) {
    int nr = Qw.GetNrows();
    Qw.ResizeTo(nr+ndiff, Qw.GetNcols()+ndiff);
    for (int j=0; j<ndiff; j++)
      Qw(nr+j, nr+j) = 1.;
  }
  if (debug) {
    cout << "Qw: ";  Qw.Print();
  }
  
  U = U*Qw;

  if (debug) {
    cout << "U: ";  U.Print();
  }
  
  // Finally, assign C and S matrices
  if (fVerbosity)
    Info("UnfoldingUtils::CSDecompQ1Taller()","Formatting output...");
  for (int j=0; j<q1; j++) 
    C(j,j) = alpha(j);
  for (int i=0; i<q2; i++) 
    S(i,l-p+i) = beta(l-p+i);

  if (debug) {
    cout << "alpha: ";  alpha.Print();
    cout << "beta: ";   beta.Print();
    cout << "C: ";  C.Print();
    cout << "S: ";  S.Print();
    TMatrixD upper(C, TMatrixD::kMultTranspose, Z);
    TMatrixD lower(S, TMatrixD::kMultTranspose, Z);
    upper = U*upper;
    upper = Q1-upper;
    lower = V*lower;
    lower = Q2-lower;
    cout << "Q1 - U*C*Z\': ";  upper.Print();
    cout << "Q2 - V*S*Z\': ";  lower.Print();
  }

  cs.C.ResizeTo(C);
  cs.S.ResizeTo(S);
  cs.U.ResizeTo(U);
  cs.V.ResizeTo(V);
  cs.Z.ResizeTo(Z);
  cs.alpha.ResizeTo(l);
  cs.beta.ResizeTo(l);

  cs.C = C;
  cs.S = S;
  cs.U = U;
  cs.V = V;
  cs.Z = Z;
  cs.alpha = alpha;
  cs.beta = beta;

  return cs;
}


GSVDecompResult 
UnfoldingUtils::GSVD(TMatrixD& A, TMatrixD& B)
{
  bool debug = false;
  GSVDecompResult g;
  int m,p,n,r;
  m = A.GetNrows();
  n = A.GetNcols();
  p = B.GetNrows();

  // TODO: check # A cols = B cols, and  m >= n >= p

  TMatrixD M(m+p,n);

  // 1. SVD of M = [A;B]: M = Q*S*Z'
  // TMatrixDSub(M, 0, m-1,   0, n-1) += A;
  // TMatrixDSub(M, m, m+p-1, 0, n-1) += B;
  M.SetSub(0, 0, A);
  M.SetSub(m, 0, B);

  if (M.GetNoElements() > 100000) 
    Printf("UnfoldingUtils::GSVD(): " 
	   "Computing initial SVD on M (%d x %d)...", m+p, n);
  TDecompSVD svdM(M);
  TMatrixD Q  = svdM.GetU(); // m+p x m+p
  TMatrixD Z  = svdM.GetV(); //   n x n
  TVectorD sv = svdM.GetSig();

  r = sv.GetNrows(); // Assume M has full rank to save time
  //r = Rank(M);

  if (debug)
    Printf("Rank(M): %d. M = QSZ\', Q (%d x %d) Z (%d x %d)", 
	   r, Q.GetNrows(), Q.GetNcols(), Z.GetNrows(), Z.GetNcols());

  TMatrixD Sr(r,r); // Sigma_r submatrix from Sigma
  for (int i=0; i<r; i++) 
    Sr(i,i) = sv(i);

  // 2. Partition Q to match dimensions of A and B.
  TMatrixD Q1 = Q.GetSub(0,m-1,0,n-1);
  TMatrixD Q2 = Q.GetSub(m, m+p-1, 0, n-1);  
  bool q1Taller = Q1.GetNrows() > Q2.GetNrows();

  // 3. Do CS decomposition
  CSDecompResult csd = 
    (q1Taller)? CSDecompQ1Taller(Q1, Q2) : CSDecomp(Q1, Q2);
  
  // 4. Assign output struct members
  g.U.ResizeTo(m,n); //csd.U);
  g.V.ResizeTo(csd.V);
  g.C.ResizeTo(n,n); //csd.C);
  g.S.ResizeTo(csd.S);
  g.alpha.ResizeTo(csd.alpha);
  g.beta.ResizeTo(csd.beta);
  g.gamma.ResizeTo(n);

  g.U = csd.U.GetSub(0,m-1,0,n-1);
  g.V = csd.V;
  g.C = csd.C.GetSub(0,n-1,0,n-1);
  g.S = csd.S;
  g.alpha = csd.alpha;
  g.beta = csd.beta;
  g.gamma = ElemDiv(g.alpha, g.beta, 9999e12);
  
  // Create X' from V'Sigma (upper left) and W (lower right)
  TMatrixD XT(n,n);
  XT.SetSub(0,0,TMatrixD(csd.Z,TMatrixD::kTransposeMult,Sr));
  for (int i=r; i<n; i++) {
    XT(i,i) = Sr(r-1,r-1);
  }
  XT = TMatrixD(XT, TMatrixD::kMultTranspose, Z);

  g.XT.ResizeTo(XT);
  g.XT = XT;

  return g;  
}


void 
UnfoldingUtils::ReverseColumns(TMatrixD& A, int col1, int col2)
{
  // Reverse the column sequence col1...col2, col1 and col2 included.
  TMatrixD B(A);
  int n = B.GetNcols();
  int ncols = col2-col1+1;
  if (ncols > n) {
    Error("ReverseColumns()",
	  "Requested column range (%d-%d) > out of bounds (0-%d).",
	  col1,col2,A.GetNcols());
    return;
  }

  for (int j=col1; j<=col2; j++) { // j is col index of B
    TMatrixDColumn(B, j) = TMatrixDColumn(A, ncols-j-1);
  }
  A = B;
}

void 
UnfoldingUtils::ReverseColumns(TMatrixD& A)
{
  TMatrixD B(A);
  int n = B.GetNcols();

  for (int j=0; j<n; j++) {
    TMatrixDColumn(B, j) = TMatrixDColumn(A, n-j-1);
  }
  A = B;
}

void 
UnfoldingUtils::ReverseRows(TMatrixD& A)
{
  TMatrixD B(A);
  int m = B.GetNrows();

  for (int i=0; i<m; i++) {
    TMatrixDColumn(B, i) = TMatrixDColumn(A, m-i-1);
  }
  A = B;
}

void 
UnfoldingUtils::ReverseVector(TVectorD& v)
{
  TVectorD vtmp(v);
  int m = v.GetNrows();
  for (int i=0; i<m; i++) {
    vtmp(i) = v(m-i-1);
  }
  v = vtmp;
}

void 
UnfoldingUtils::SwapColumns(TMatrixD &A, int col1, int col2)
{
  int nc = A.GetNcols();
  if (col1 >= nc || col2 >= nc)
    Error("SwapColumns", "col1 or col2 index out of bounds");

  TMatrixD B(A);
  
  TMatrixDColumn(B, col1) = TMatrixDColumn(A, col2);
  TMatrixDColumn(B, col2) = TMatrixDColumn(A, col1);

  A = B;
}

void 
UnfoldingUtils::SwapElements(TVectorD& v, int j1, int j2)
{
  int nr = v.GetNrows();
  if (j1 >= nr || j2 >= nr)
    Error("SwapElements", "an index is out of bounds");
  
  TVectorD v2(v);
  
  v2(j1) = v(j2);
  v2(j2) = v(j1);
  
  v = v2;
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
