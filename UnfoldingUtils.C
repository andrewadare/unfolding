#include "UnfoldingUtils.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
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
  fHistEff(0)
{
}

UnfoldingUtils::UnfoldingUtils(TH2D* hA, 
			       TH1D* hMeas, 
			       TH2D* hMeasCov, 
			       TH1D* hXini, 
			       TH1D* hXtrue, 
			       TH1D* hEff) :
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
  fHistEff(hEff)
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
  fVecbTilde.ResizeTo(fM);
  fVecXini.ResizeTo(fN);
  fVecXtrue.ResizeTo(fN);

  fMatA = Hist2Matrix(fHistA);
  fVecb = Hist2Vec(fHistMeas);

  if (fHistXini)
    fVecXini = Hist2Vec(fHistXini);
  else
    for (int j=0; j<fN; j++)
      fVecXini(j) = 1.0;
  
  if (fHistXtrue)
    fVecXtrue = Hist2Vec(fHistXtrue);

  // Probability matrix Ahat
  if (!fHistAProb) {
    fHistAProb = (TH2D*) fHistA->Clone("fHistAProb");
    NormalizeXSum(fHistAProb, (fHistEff) ? fHistEff : 0);
  }
  fMatAhat = Hist2Matrix(fHistAProb);

  // Data uncertainty
  TVectorD eb(fM);  
  for (int i=0; i<fM; i++)
    eb(i) = fHistMeas->GetBinError(i+1);
  
  // Data covariance matrix
  if (fHistMeasCov)
    fMatB = Hist2Matrix(fHistMeasCov);
  else
    for (int i=0; i<fM; i++) {
      fMatB(i,i) = eb(i)*eb(i);
    }
  fMatBinv = MoorePenroseInverse(fMatB);

  // Create rescaled (~) quantities
  if (0) {
    // Hocker eq. (33)
    // Decompose B = QRQ' where R_{ij} = r_i^2 \delta_{ij}
    TDecompSVD QRQT(fMatB); 
    TMatrixD Q = QRQT.GetU(); // Q=U=V if B is symm. & pos.-semidef
    TVectorD rsq = QRQT.GetSig(); // R(i,i) \equiv rsq(i)
    
    // Hocker eq. (34)
    fMatATilde = Q*fMatA;
    for (int i=0;i<fM;i++) // row index 
      for (int j=0;j<fN;j++) // col index
	fMatATilde(i,j) *= 1./TMath::Sqrt(rsq(i));
    fVecbTilde = Q*fVecb;
    for (int i=0;i<fM;i++)
      fVecbTilde(i) *= 1./TMath::Sqrt(rsq(i));
  }
  
  // For now, assume uncorrelated errors in b

  fMatATilde = DivColsByVector(fMatA, eb);
  fHistATilde = Matrix2Hist(fMatATilde, "fHistATilde",
			    fMeasX1,fMeasX2,fTrueX1,fTrueX2);
  fVecbTilde = ElemDiv(fVecb, eb);
  fHistbTilde = Vec2Hist(fVecbTilde, fMeasX1, fMeasX2, 
			 "fHistbTilde", "Scaled measured distribution");
  for (int i=0;i<fM;i++)
    fHistbTilde->SetBinError(i+1, 1.0);
  return;
}

bool
UnfoldingUtils::BinningOk()
{
  // Check for binning incompatibilities
  bool isok = true;
  
  int nMeas = fHistMeas->GetNbinsX(); // Must equal fM
  if (nMeas != fM) {
    gROOT->Warning("", "Meas. bin mismatch: hMeas %d, TH2 (x-axis) %d",
		   nMeas, fM);
    isok = false;
  }
  if (fHistXini) {
    int nXini = fHistXini->GetNbinsX(); // Must equal fN
    if (nXini != fN) {
      gROOT->Warning("", "True bin mismatch: x^ini %d, TH2 (y-axis) %d",
		     nXini, fN);
      isok = false;
    }
  }
  if (fHistXtrue) {
    int nXtrue = fHistXtrue->GetNbinsX(); // Must equal fN
    if (nXtrue != fN) {
      gROOT->Warning("", "True bin mismatch: hXtrue %d, TH2 (y-axis) %d",
		     nXtrue, fN);
      isok = false;
    }
  }
  return isok;
}

void 
UnfoldingUtils::SetTrueRange(double x1, double x2)
{
  fTrueX1 = x1; 
  fTrueX2 = x2;
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
  TMatrixD A(fMatA);
  TVectorD b(fVecb);
  if (fTilde) {
    A = fMatATilde;
    b = fVecbTilde;
  }
  else if (opt.Contains("^")) {
    A = fMatAhat;
  }

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

GSVDResult 
UnfoldingUtils::GSVDAnalysis(TMatrixD& L, double lambda, TH2* hA, TH1* hb, TString opt)
{
  // Decompose A, L jointly as A = U*C*X', L = V*S*X'. 
  // If A is m x n, and L is p x n, and A has full rank,
  // U  is m x m 
  // V  is p x p
  // X' is n x n
  // C  is m x n
  // S  is p x n
  // alpha and beta are the "interesting" diagonal elements of C,S
  // respectively, and have length p.

  static int id = 0; id++;
  GSVDResult result;

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

  TVectorD regc = ElemMult(f,c);
  TVectorD wreg = X*regc;
  TVectorD xreg = ElemMult(fVecXini, wreg);

  // Save to output for parameter optimization analysis
  result.bInc.ResizeTo(m);
  result.bInc = (LMatrix(m,kUnitMatrix) - g.U * UT)*b;

  // Assign output struct members
  result.n      = n;
  result.m      = m;
  result.p      = p;
  result.lambda = lambda;
  result.alpha.ResizeTo(n);  result.alpha  = g.alpha;
  result.beta .ResizeTo(n);  result.beta   = g.beta;
  result.gamma.ResizeTo(n);  result.gamma  = g.gamma;
  result.f    .ResizeTo(n);  result.f      = f;
  result.UTb  .ResizeTo(n);  result.UTb    = utb;
  result.coeff.ResizeTo(n);  result.coeff  = c;
  result.regc .ResizeTo(n);  result.regc   = regc;

  // Copy A,L,and b to output to ensure that the correct quantities
  // are used later.
  result.X.ResizeTo(X);      result.X = X;
  result.U.ResizeTo(g.U);    result.U = g.U;
  result.V.ResizeTo(g.V);    result.V = g.V;
  result.L.ResizeTo(L);      result.L = L;
  result.A.ResizeTo(A);      result.A = A;
  result.b.ResizeTo(b);      result.b = b;

  result.UHist  = Matrix2Hist(g.U, Form("U_gsvd_%d",id),
			      fTrueX1, fTrueX2,0,m);
  result.XHist  = Matrix2Hist(X, Form("X_gsvd_%d",id),
			      fTrueX1, fTrueX2,0,m);
  result.wregHist = Vec2Hist(wreg, fTrueX1,fTrueX2,
			     Form("gsvd_wreg_%d",id), 
			     Form("w (#lambda = %g)", lambda));
  result.xregHist = Vec2Hist(xreg, fTrueX1,fTrueX2,
			     Form("gsvd_xreg_%d",id),
  			     Form("x (#lambda = %g)", lambda));
  // Absolute values for plotting
  TVectorD utbAbs(utb);
  TVectorD cAbs(c);
  TVectorD rcAbs(regc);
  for (int i=0; i<n; i++) {
    if (utbAbs(i) < 0) utbAbs(i) *= -1;
    if (cAbs(i) < 0)     cAbs(i) *= -1;
    if (rcAbs(i) < 0)   rcAbs(i) *= -1;
  }
  result.UTbAbs   = Vec2Hist(utbAbs, 0,n,Form("gsvd_utb_abs%d",id), "#||{u^{T}#upointb} ");
  result.coeffAbs = Vec2Hist(cAbs, 0,n,Form("gsvd_c_abs%d",id), "#||{u^{T}#upointb}/#alpha ");
  result.regcAbs  = Vec2Hist(rcAbs, 0,n,Form("gsvd_rc_abs%d",id), "f#||{u^{T}#upointb}/#alpha ");

  SetTH1Props(result.UTbAbs,   kBlue, 0, kBlue, kFullSquare, 1.0);
  SetTH1Props(result.coeffAbs, kRed, 0, kRed, kOpenSquare, 1.0);
  SetTH1Props(result.regcAbs, kMagenta+1, 0, kMagenta+1, kOpenCircle, 1.0);
  SetTH1Props(result.xregHist, kGreen+2, 0, kGreen+2, kFullCircle, 1.0);

  return result;
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
UnfoldingUtils::DrawGSVDPlot(TObjArray* svdhists, double ymin, double ymax, TString opt)
{
  static int i=0; i++;
  TCanvas* c = new TCanvas(Form("cgsvd%d",i), Form("cgsvd%d",i), 1);
  TH1D* hs = (TH1D*)svdhists->At(0); // Sing. values
  TH1D* hd = (TH1D*)svdhists->At(1); // d (lambda=0)
  TH1D* hl = (TH1D*)svdhists->At(2); // d (lambda)
  if (!hs) Warning("DrawGSVDPlot()","!hs");
  if (!hd) Warning("DrawGSVDPlot()","!hd");
  if (!hl) Warning("DrawGSVDPlot()","!hl");
  if (hs) { // Draw frame histo
    TH2F* hsvd = new TH2F(Form("hgsvd%d",i), "GSVD components;column index i;", 
			  200, 0, hs->GetNbinsX(), 200, ymin, ymax);
    TLine l;
    l.SetLineColor(kGray);

    hsvd->Draw();
    l.DrawLine(0, 1, (double)hs->GetNbinsX(), 1);
    if (opt.Contains("hs"))
      hs->Draw("plsame");
    hd->Draw("plsame");
    hl->Draw("plsame");
    gPad->SetLogy();
    
    TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
    if (opt.Contains("hs"))
      leg->AddEntry(hs, hs->GetTitle(), "p");
    leg->AddEntry(hd, hd->GetTitle(), "ep");
    leg->AddEntry(hl, hl->GetTitle(), "ep");
    leg->SetFillColor(kNone);
    leg->Draw();
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
  
  TestProblem t;
  double dm = (xm2-xm1)/m;
  double dt = (xt2-xt1)/n;
  
  // Create response matrix from kernelFn, discretized via simple colocation. 
  TH1D* hKernel = new TH1D("hKernel", "Discretized convolution kernel", 
			   m+n+1, -xt2-dt/4, xm2+dm/4);
  for (int j=1; j<=2*m+1; j++) {
    double ctr = hKernel->GetBinCenter(j);
    double val = kernelFn->Eval(ctr);
    hKernel->SetBinContent(j, val);
  }
  t.Response = BandedDiagonalMatrix(hKernel, m, n, xt1, xt2, xm1, xm2);
  t.Response->SetTitle(Form("%d x %d convolution matrix;"
			    "s (observed);t (true)", m, n));
  
  // Model a true and a measured distribution
  // There is no bIdeal for this problem.
  t.xTruth    = new TH1D("hTrue",     "",     n, xt1, xt2);  
  t.xTruthEst = new TH1D("hTrueEst", "hTrueEst", n, xt1, xt2);  
  t.bNoisy    = new TH1D("hMeas",     "hMeas",     m, xm1, xm2);  
  for (Int_t i=0; i<nEvents; i++) {
    Double_t xt = truthFn->GetRandom();
    t.xTruthEst->Fill(xt);
    Double_t xsmear = kernelFn->GetRandom();
    t.bNoisy->Fill(xt + xsmear);
  }

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
UnfoldingUtils::Hist2Vec(const TH1* h)
{
  // Returns TVectorD of the bin contents of the input histogram
  int nb = h->GetNbinsX();
  TVectorD v(nb);
  if (!h) return v;
  for (Int_t i= 0; i<nb; i++) {
    v(i) = h->GetBinContent(i+1);
  }
  return v;
}

double
UnfoldingUtils::SmoothingNorm(TVectorD& x, int regtype)
{
  double sn = 0;
  switch (regtype) {
  case kNone:     sn = 0; break;
  case kTotCurv:  sn = Curvature(x); break;
  }
  return sn;
}

double
UnfoldingUtils::Curvature(TVectorD& x)
{
  // Eq. (38), NIM A 372 (1996) 469-481. 
  double delta=0, val=0;
  for (int i=1; i<fN-1; i++) {
    delta = (x(i+1) - x(i)) - (x(i) - x(i-1));
    val += delta*delta;
  }

  return val;
}

double
UnfoldingUtils::RegChi2(const double *pars)
{
  double chi2 = 0;
  // Fit parameters vector
  TVectorD x(fN);
  for (int i=0; i<fN; i++)
    x(i) = pars[i];
  
  // Unmodified chi^2 (Ax-b)'*Binv*(Ax-b) (Hocker eq. (30))
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

TH1D*
UnfoldingUtils::UnfoldChiSqMin(double regWt,
			       TObjArray* /*output*/,
			       TString opt,
			       TH1* hXStart,
			       TH2* hA,
			       TH1* hb,
			       TH1* hXini)
{
  static int id=0; id++;
  
  fRegWeight = regWt;
  fTilde = (opt.Contains("~")) ? true : false;

  Info("UnfoldingUtils::UnfoldChiSqMin()",
       "Using reg type %d, weight %g. Uncertainty rescaling: %s",
       fRegType, fRegWeight, fTilde? "yes":"no");
  
  if (!hA)
    hA = (fTilde) ? fHistATilde : fHistAProb;
  if (!hb)
    hb = (fTilde) ? fHistbTilde : fHistMeas;
  if (fTilde && !hXini)
    hXini = fHistXini;
  
  // Set up unfolded result. 
  // TH1 bin info from truth axis of response matrix.
  int nBinsT = hA->GetNbinsY();           // # true/gen/unfolded bins
  int nPars = nBinsT;                     // Add 1 for not-found events (?)
  double xt1 = hA->GetYaxis()->GetXmin();
  double xt2 = hA->GetYaxis()->GetXmax();
  TH1D* hUnf = new TH1D(Form("hChsq%d",id), Form("hChsq%d",id),
			nBinsT, xt1, xt2);

  // Starting values for the minimizer. If not passed in, use hXini. If
  // no hXini, use measured points.
  if (!hXStart) {
    if (hXini)
      hXStart = (TH1D*)hXini->Clone(Form("hTMinStart%d",id));
    else {
      hXStart = (TH1D*)hUnf->Clone(Form("hTMinStart%d",id));
      for (int i=0; i<nBinsT; i++) {
	int bin = hb->FindBin(hXStart->GetBinCenter(i+1));
	hXStart->SetBinContent(i+1, hb->GetBinContent(bin));
      }
    }
  }

  // Set up the chi^2 minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.0001);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(this, &UnfoldingUtils::RegChi2, nPars);
  min->SetFunction(f);
  
  // Initialize tmx array and pass to minimizer.
  double* tmx = new double[nBinsT];
  double stepSize = 0.1;

  for (int j=0; j<nBinsT; j++) {
    tmx[j] = hXStart->GetBinContent(j+1);

    if (hXini) // scale to w
      tmx[j] /= hXini->GetBinContent(j+1);
    
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

  Info("UnfoldingUtils::UnfoldChiSqMin()", 
       "Initial (regularized) chi squared = %g", RegChi2(tmx));
  min->Minimize(); 
  
  for (int j=0; j<nPars; j++) {
    tmx[j] = min->X()[j];
    double val = tmx[j];
    double err = val < 1 ? 1. : min->Errors()[j];
    hUnf->SetBinContent(j+1, val);
    hUnf->SetBinError(j+1, err);
  }

  TVectorD x = Hist2Vec(hUnf);
  
  // Scale w --> x. w = x if "~" option is not set.
  if (hXini)
    hUnf->Multiply(hXini);
  
  Info("UnfoldingUtils::UnfoldChiSqMin()", 
       "Final (regularized) chi squared = %g, curvature = %g",
       RegChi2(tmx), Curvature(x));
  return hUnf;
}

UnfoldingResult
UnfoldingUtils::UnfoldTikhonovGSVD(GSVDResult& gsvd,
				   TVectorD& lambda, 
				   TString /*opt*/)
{
  UnfoldingResult result;
  static int id=0; id++; // So this fn. can be called more than once

  int nk = lambda.GetNrows();
  int n  = gsvd.n;
  int m  = gsvd.m; 
  int p  = gsvd.p;
  result.XReg.ResizeTo(n, nk);
  result.LCurve = new TGraph(nk);
  result.GcvCurve = new TGraph(nk);
  result.lambdaGcv = 0;
  double gcvMin = 1e99;

  // Scan over lambda values, generate nk solutions
  for (int k=0; k<nk; k++) {
    
    // Create Tikhonov filter factors for this lambda
    double l = lambda(k);
    TVectorD f(n);
    for (int i=0; i<n; i++) {
      if (i >= n-p) {
	double g2 = gsvd.gamma(i)*gsvd.gamma(i); 
	f(i) = g2 / (g2 + l*l);
      }
      else
	f(i) = 1.0;
    }
    
    // Damped GSVD coefficients
    TVectorD regc = ElemMult(f, gsvd.coeff);
    TVectorD wreg = gsvd.X*regc;
    TVectorD xreg = ElemMult(fVecXini, wreg);

    // Regularized solution
    for (int j=0; j<n; j++)
      result.XReg(j,k) = xreg(j);

    // Parameter optimization analysis -------------------------------
    // ---------------------------------------------------------------

    // Compute Lx_reg
    TVectorD vec = ElemDiv(gsvd.UTb, gsvd.gamma);
    vec = ElemMult(f, vec);
    TVectorD Lx = vec.GetSub(n-p, n-1);
    Lx = gsvd.V * Lx;
    double lxnorm = TMath::Sqrt(Lx*Lx);

    // Compute r = b - Ax_reg
    TVectorD f1(f);
    for (int j=0; j<n; j++)
      f1(j) = 1-f(j);
    vec = ElemMult(f1, gsvd.UTb);
    TMatrixD Up = gsvd.U.GetSub(0,m-1,m-p,m-1); // m x p
    TVectorD r = Up * vec.GetSub(n-p, n-1) - gsvd.bInc;
    double rnorm = TMath::Sqrt(r*r);

    double gcv = rnorm / (m - f.Sum());

    result.LCurve->SetPoint(k, rnorm, lxnorm);
    result.GcvCurve->SetPoint(k, lambda(k), gcv);
    if (gcv < gcvMin) {
      gcvMin = gcv;
      result.lambdaGcv = lambda(k);
      result.kGcv = k;
    }

  }


  result.LCurve->SetTitle("GSVD L-Curve;||Ax_{#lambda}-b||_{2};||Lx_{#lambda}||_{2}");
  result.GcvCurve->SetTitle("GSVD cross-validation curve;"
			    "#lambda;G(#lambda)");

  result.XRegHist = Matrix2Hist(result.XReg, Form("xregHist_%d",id),
				fTrueX1, fTrueX2,lambda(0),lambda(nk-1));
  result.XRegHist->SetTitle("GSVD solutions;x;#lambda");
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
UnfoldingUtils::UnfoldRichardsonLucy(const int nIterations, 
				     TString opt, 
				     const TH1* hXStart)
{
  // See eq. (4.3), J. Bardsley & C. Vogel, 
  // SIAM J. SCI. COMPUT. Vol. 25, No. 4, pp. 1326–1343
  
  UnfoldingResult result;
  static int ihist = 0; ihist++;       // Unique ID

  result.LCurve = new TGraph(nIterations);
  result.GcvCurve = new TGraph(nIterations);
  result.lambdaGcv = 0;
  //  double gcvMin = 1e99;
  result.LCurve->SetNameTitle(Form("LCurve_RL_%d",ihist), 
			      Form("Richardson-Lucy L-Curve;"
				   "Residual norm ||Ax_{k}-b||_{2};"
				   "Solution norm ||x_{k}||_{2}"));
  TMatrixD A = GetA(opt);
  TVectorD b = Getb(opt);
  TVectorD xini(fVecXini);
  TVectorD ones(fM); 
  for (int i=0; i<fM; i++) 
    ones(i)=1.0;

  result.WReg.ResizeTo(fN, nIterations);
  result.XReg.ResizeTo(fN, nIterations);

  // Prior vector x0 must be positive
  TVectorD x0 = (hXStart)? Hist2Vec(hXStart) : ones;
  for (int j=0; j<x0.GetNrows(); j++) {
    if (x0(j)<=0) {
      Warning("UnfoldingUtils::UnfoldRichardsonLucy()",
	      "Initial point x0(%d) = %g must be positive, "
	      "setting to 1.0", j, x0(j) );
      x0(j) = 1.;
    }
  }
  
  double x1 = fTrueX1, x2 = fTrueX2;
  double hx1 = hXStart->GetXaxis()->GetXmin();
  double hx2 = hXStart->GetXaxis()->GetXmax();
  if (hx1 != x1)
    Warning("UnfoldingUtils::UnfoldRichardsonLucy()",
	    "hXStart x1 %g != stored x1 %g", hx1, x1);
  if (hx2 != x2)
    Warning("UnfoldingUtils::UnfoldRichardsonLucy()",
	    "hXStart x2 %g != stored x2 %g", hx2, x2);
  
  TMatrixD AT(TMatrixD::kTransposed, A);
  TVectorD x(x0);
  TVectorD bvar(fM); // Variance of meas. data
  for (int i=0; i<fM; i++)
    bvar(i) = fMatB(i,i);
  
  b += bvar;

  for (int k=0; k<nIterations; k++) {
    printf("Richardson-Lucy iteration %d\r", k+1);
  
      TVectorD Ax = A*x + bvar; 
      TVectorD r = ElemDiv(b, Ax);
      TVectorD ATr = AT*r;
      TVectorD AT1 = AT*ones; // efficiency correction factor
      TVectorD xoverAT1 = ElemDiv(x, AT1);
      x = ElemMult(xoverAT1, ATr);
      TVectorD resid = A*x-b;
      double rnorm = TMath::Sqrt(resid*resid);
      double xnorm = TMath::Sqrt(x*x);
  
      for (int j=0; j<fN; j++)
	result.WReg(j,k) = x(j);
      
      if (fTilde)
	x = ElemMult(x,xini);
      
      for (int j=0; j<fN; j++)
	result.XReg(j,k) = x(j);
      
      result.LCurve->SetPoint(k,rnorm,xnorm);
  } // end iteration loop
  cout << endl;

  result.XRegHist = Matrix2Hist(result.XReg, Form("X_RL_%d",ihist),
				fTrueX1,fTrueX2,0,nIterations);
  
  /*
  result.wregHist = Vec2Hist(wreg, fTrueX1,fTrueX2,
  			     Form("gsvd_wreg_%d",id), 
  			     Form("w (#lambda = %g)", lambda));
  result.xregHist = Vec2Hist(xreg, fTrueX1,fTrueX2,
			     Form("gsvd_xreg_%d",id),
			     Form("x (#lambda = %g)", lambda));
  */

  // TH1D* hx=0;
  // if (hists)
  //   hx = (TH1D*)hists->At(nIterations);
  // else {
  //     hx = XHist(x,Form("RL%d",ihist),nIterations,x1,x2,0.,"");
  // }
  // if (extras) {
  //   if (gL) extras->Add(gL);
  // }

  return result;
}

TH1D*
UnfoldingUtils::UnfoldPCGLS(const int nIterations, 
			    TObjArray* hists,
			    TObjArray* extras,
			    int LType,
			    TString opt,
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

  static int ihist = 0; ihist++;       // Unique ID
  fTilde = (opt.Contains("~")) ? true : false;

  TMatrixD A(fTilde? fMatATilde : fMatA);
  TVectorD b(fTilde? fVecbTilde : fVecb);
  if (hA) A = Hist2Matrix(hA);
  if (hb) b = Hist2Vec(hb);

  TMatrixD L = LMatrix(A.GetNcols(), LType, 1e-5);
  TMatrixD W = Null(L);          // Nullspace of L

  int ncols = A.GetNcols();
  int p = L.GetNrows();
  double x1 = fTrueX1, x2 = fTrueX2;
  TVectorD xini(fVecXini);
  if (hXini) {
    x1 = hXini->GetXaxis()->GetXmin();
    x2 = hXini->GetXaxis()->GetXmax();
    xini = Hist2Vec(hXini);
  }
  
  TGraph* gL = 0;
  if (extras) {
    gL = new TGraph(nIterations);
    gL->SetLineColor(kRed);
    gL->SetMarkerColor(kRed);
    gL->SetMarkerStyle(kFullCircle);
    gL->SetMarkerSize(1.5);
    gL->SetLineWidth(2);
    gL->SetNameTitle(Form("LCurve_PCGLS_%d",ihist), 
		     Form("PCGLS L-Curve;"
			  "Residual norm ||Ax_{k}-b||_{2};"
			  "Solution norm ||Lx_{k}||_{2}"));
  }

  if (0) {  
    Printf("A (%d x %d)", A.GetNrows(), ncols);
    Printf("L (%d x %d)", p, L.GetNcols());
    Printf("W (%d x %d)", W.GetNrows(), W.GetNcols());
  }
  
  // Store q1 vectors as columns of Q1n for re-orthogonalization
  TMatrixD Q1n(p, nIterations+1);

  // Prepare for computations with L_p (= pinit.m):
  // T <-- pinv(AW)*A; x0 <-- W*pinv(AW)*b
  TMatrixD T(0,0);
  TVectorD x0(ncols); x0.Zero();

  // Dimension of nullspace of L
  int nullity = W.GetNcols();

  if (nullity>0) {
    TMatrixD AW(A*W);
    if (1) {
      Printf("AW (%d x %d)", AW.GetNrows(), AW.GetNcols());    
    }
    TMatrixD AWinv = MoorePenroseInverse(AW);
    if (0) {
      Printf("MoorePenroseInverse(A*W) (%d x %d)",
	     AWinv.GetNrows(), AWinv.GetNcols());
      AWinv.Print();
    }
    T.ResizeTo(AWinv.GetNrows(), ncols);
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
  
  for (int n=0; n<=nIterations; n++) {
    printf("PCGLS iteration %d\r", n);
    
    // Run CGLS algorithm for n=1,2,...
    if (n>0) {
      TVectorD Az = A*z;
      alpha = dq / Az.Norm2Sqr();
      x += alpha*z;
      r -= alpha*Az;
      s  = AT*r;
      LTSolve(q1, L, s);

      // Reorthogonalize q1 to previous q1 vectors
      for (int i=0; i<n; i++) {
	TVectorD qi = TMatrixDColumn(Q1n, i);
	q1 -= (qi*q1)*qi; // Modified Graham-Schmidt
	q1norm = TMath::Sqrt(q1*q1);
	if (q1norm > 0) q1norm = 1./q1norm;
	TMatrixDColumn(Q1n, n) = q1norm*q1;
      }
      
      LSolve(q, L, q1, W, T);
      dq2 = s*q;
      beta = dq2/dq;
      dq = dq2;
      z = q + beta*z;

      if (gL)    // Set L-Curve points
	gL->SetPoint(n-1, TMath::Sqrt(r*r), TMath::Sqrt((L*x)*(L*x)));
    }
    
    if (hists) { // Save nth iterate as a TH1D
      if (hXini) {
	TVectorD xrescaled(ncols);
	for (int j=0;j<ncols;j++) xrescaled(j) = x(j)*xini(j);
	hists->Add(XHist(xrescaled,Form("PCGLS%d",ihist),n,x1,x2,0.,""));
      }
      else
	hists->Add(XHist(x,Form("PCGLS%d",ihist),n,x1,x2,0.,""));
    }
  } // end iteration loop
  cout << endl;
  
  TH1D* hx=0;
  if (hists)
    hx = (TH1D*)hists->At(nIterations);
  else {
    if (hXini) {
      TVectorD xrescaled(ncols);
      for (int j=0;j<ncols;j++) xrescaled(j) = x(j)*xini(j);
      hx = XHist(xrescaled,Form("PCGLS%d",ihist),nIterations,x1,x2,0.,"");
    }
    else
      hx = XHist(x,Form("PCGLS%d",ihist),nIterations,x1,x2,0.,"");
  }
  if (extras) {
    if (gL) extras->Add(gL);
  }

  return hx;
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
  // base should have a form like "CGLS"+"iHist"
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
    Error("BkgMatrix()", "nNegBins < 1 (%d)", nNegBins);
  if (nMeas > nPosBins)
    Error("BkgMatrix()", "nMeas (%d) > nPosBins (%d)", nMeas, nPosBins);
  if (nTrue >= zerobin)
    Error("BkgMatrix()", "nTrue (%d) > nNegBins (%d)", nTrue, nNegBins);
  
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
				TString opt)
{
  // Propagate the measured data cov matrix to the unfolded result.
  static int id=0; id++;
  int seed = 0; // 0 means use TUUID identifier
  TRandom3 rand(seed);
 
  // Measurement covariance 
  TMatrixD B(fMatB);
  if (fTilde)
    B = LMatrix(fM, kUnitMatrix); // B~ is unit matrix
  
  TDecompChol decomp(fMatB);
  TMatrixD UT(TMatrixD::kTransposed, decomp.GetU());
  TH1D* hUnf = 0; // Unfolded result

  // Returned covariance matrix
  TString title("Covariance matrix from MC trials");
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
    for (int i=0; i<fM; i++) {
      g(i) = rand.Gaus(0,1);
    }
    
    g *= B;//UT;
    for (int i=0; i<fM; i++) {
      hbTrial->SetBinContent(i+1, fVecb(i) + g(i));
      hbTrial->SetBinError(i+1, fHistMeas->GetBinError(i+1));
    }
    
    switch (algo) 
      {
      case kSVDAlgo: 
	hUnf = UnfoldSVD(regPar, 0, opt, 0, hbTrial, 0); 
	break;
	// case kRichLucyAlgo: 
	// 	h = UnfoldRichardsonLucy(); 
	// 	break;
	// case kChi2MinAlgo: 
	// 	h = UnfoldChi2Min(); 
	// 	break;
	// case kPCGLSAlgo: 
	// 	h = UnfoldPCGLS(); 
	// 	break;
      }
    
    for (int j=1; j<=fN; j++)
      hxMean->AddBinContent(j,hUnf->GetBinContent(j)/nTrials); 
  }

  rand.SetSeed(seed);
  
  // Loop 2: Compute E[(x-E[x])*(x-E[x])'] as hcov
  for (int t=0; t<nTrials; t++) {  
    for (int i=0; i<fM; i++) {
      g(i) = rand.Gaus(0,1);
    }
    g *= B;//UT;
    for (int i=0; i<fM; i++) {
      hbTrial->SetBinContent(i+1, fVecb(i) + g(i));
      hbTrial->SetBinError(i+1, fHistMeas->GetBinError(i+1));
    }

    switch (algo) 
      {
      case kSVDAlgo: 
	hUnf = UnfoldSVD(regPar, 0, opt, 0, hbTrial, 0); 
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
  TDecompSVD svdQ1(Q1);
  TMatrixD U     = svdQ1.GetU();    // m x m
  TMatrixD Z     = svdQ1.GetV();    // l x l
  alpha          = svdQ1.GetSig();  // l

  // 3.
  for (int i = 0; i<l; i++) {
    beta(i)  = 1.; 
    if (i>=q1) alpha(i) = 0.;
  }
  
  Printf("CSD (Q1 taller): m=%d, p=%d, l=%d, q1=%d, q2=%d",  m,p,l,q1,q2);
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
  for (int j=0; j<q1; j++) 
    C(j,j) = alpha(j);
  for (int i=0; i<q2; i++) 
    S(i,l-p+i) = beta(l-p+i);

    cout << "alpha: ";  alpha.Print();
    cout << "beta: ";   beta.Print();
  
  if (debug) {
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
  TMatrixDSub(M, 0, m-1,   0, n-1) += A;
  TMatrixDSub(M, m, m+p-1, 0, n-1) += B;
  r = Rank(M);

  TDecompSVD svdM(M);
  TMatrixD Q  = svdM.GetU(); // m+p x m+p
  TMatrixD Z  = svdM.GetV(); //   n x n
  TVectorD sv = svdM.GetSig();

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
  g.U.ResizeTo(csd.U);
  g.V.ResizeTo(csd.V);
  g.C.ResizeTo(csd.C);
  g.S.ResizeTo(csd.S);
  g.alpha.ResizeTo(csd.alpha);
  g.beta.ResizeTo(csd.beta);
  g.gamma.ResizeTo(n);

  g.U = csd.U;
  g.V = csd.V;
  g.C = csd.C;
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
