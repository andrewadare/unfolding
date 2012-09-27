#include "UnfoldingUtils.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TDecompSVD.h"
#include "TDecompQRH.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#endif

using std::cout;
using std::endl;

ClassImp(UnfoldingUtils);

UnfoldingUtils::UnfoldingUtils() : 
  fM(0),
  fN(0),
  fMeasX1(0),
  fMeasX2(0),
  fTrueX1(0),
  fTrueX2(0),
  fRegWeight(0.),
  fRegType(kNoReg)
{
}

UnfoldingUtils::UnfoldingUtils(TH2* hAhat, TH1* hMeas, TH1* hXini, TH1* hXtrue) : 
  fM(0),
  fN(0),
  fMeasX1(0),
  fMeasX2(0),
  fTrueX1(0),
  fTrueX2(0),
  fRegWeight(0.0),
  fRegType(kNoReg)
{
  Setup(hAhat, hMeas, hXini, hXtrue);
}

int
UnfoldingUtils::Setup(TMatrixD& Ahat, TVectorD& b, TVectorD& xtrue, 
		      double mx1, double mx2, double tx1, double tx2)
{
  int status = 0; // If nonzero, problems occurred.
  fM = Ahat.GetNrows();
  fN = Ahat.GetNcols();
  fMatAhat.ResizeTo(Ahat);
  fMatAhat = Ahat;
  fVecb.ResizeTo(fM);
  fVecb = b;
  fVecxTrue.ResizeTo(fN);
  fVecxTrue = xtrue;

  fTrueX1 = tx1;
  fTrueX2 = tx2;
  fMeasX1 = mx1;
  fMeasX2 = mx2;

  // Check for binning incompatibilities
  int nMeas = b.GetNoElements();
  if (nMeas != fM) {
    gROOT->Warning("", "Meas bin mismatch: b %d, m %d", nMeas, fM);
    status++;
  }
  int nTrue = xtrue.GetNoElements();
  if (nTrue != fN) {
    gROOT->Warning("", "True bin mismatch: xtrue %d, n %d", nTrue, fN);
    status++;
  }
  
  return status;
}

int 
UnfoldingUtils::Setup(TH2* hAhat, TH1* hMeas, TH1* hXini, TH1* hXtrue)
{
  int status = 0; // If nonzero, problems occurred.
  fHistAProb = (TH2D*)hAhat->Clone("fHistAProb");
  fHistMeas = (TH1D*)hMeas->Clone("fHistMeas");
  
  fM = fHistAProb->GetNbinsX();    // Measured bins (rows)
  fN = fHistAProb->GetNbinsY();    // True/gen bins (cols)
  fTrueX1 = fHistAProb->GetYaxis()->GetXmin();
  fTrueX2 = fHistAProb->GetYaxis()->GetXmax();
  fMeasX1 = fHistAProb->GetXaxis()->GetXmin();
  fMeasX2 = fHistAProb->GetXaxis()->GetXmax();

  fMatA.ResizeTo(fM,fN);
  fMatAhat.ResizeTo(fM,fN);
  fMatATilde.ResizeTo(fM,fN);
  fMatB.ResizeTo(fM,fM);

  fVecb.ResizeTo(fM);
  fVecbTilde.ResizeTo(fM);

  // Check for binning incompatibilities
  int nMeas = hMeas->GetNbinsX(); // Must equal fM
  if (nMeas != fM) {
    gROOT->Warning("", "Meas. bin mismatch: hMeas %d, TH2 (x-axis) %d", nMeas, fM);
    status++;
  }
  if (hXini) {
    fHistXini = (TH1D*)hXini->Clone("fHistXini");
    int nXini = hXini->GetNbinsX(); // Must equal fN
    if (nXini != fN) {
      gROOT->Warning("", "True bin mismatch: x^ini %d, TH2 (y-axis) %d", nXini, fN);
      status++;
    }
  }
  if (hXtrue) {
    fHistXtrue = (TH1D*)hXtrue->Clone("fHistXtrue");
    int nXtrue = hXtrue->GetNbinsX(); // Must equal fN
    if (nXtrue != fN) {
      gROOT->Warning("", "True bin mismatch: hXtrue %d, TH2 (y-axis) %d", nXtrue, fN);
      status++;
    }
  }
  
  // "Probability" Response matrix (columns sum to 1.0)
  fMatAhat = Hist2Matrix(fHistAProb);
  
  // Copy measured hist data points to b
  for (int i=0;i<fM;i++) 
    fVecb(i)=fHistMeas->GetBinContent(i+1);
  
  // Compute B = cov(b) assuming only indep. stat. err. at this time
  // TODO: change Setup to optionally pass in a covariance
  // matrix. Then if one is not passed in, it can be computed here by
  // default.
  fMatB.UnitMatrix();
  for (int i=0;i<fM;i++) 
    fMatB(i,i) *= fHistMeas->GetBinError(i+1)*fHistMeas->GetBinError(i+1);

  fHistMeasCov = Matrix2Hist(fMatB, "fHistMeasCov",fMeasX1,fMeasX2,fMeasX1,fMeasX2);

  // Compute A, Atilde, and btilde if xini was provided
  if (hXini) {
    ComputeRescaledSystem();
  }

  if (status)
    gROOT->Warning("", "%d problem%s in UnfoldingUtils::Setup()", 
		   status, status>1 ? "s" : "");
  
  return status;
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

void 
UnfoldingUtils::ComputeRescaledSystem()
{
  // Create A by scaling probability matrix Ahat by xini, then
  // columns sum to counts (instead of 1.0). This "cancels" the
  // scaling w_j = x_j/xini_j; see A. Hocker, NIM A 372 (1996)
  // 469-481.
  for (int i=0;i<fM;i++) // row index 
    for (int j=0;j<fN;j++) // col index
      fMatA(i,j) = fMatAhat(i,j)*fHistXini->GetBinContent(j+1);
 
  fHistA = Matrix2Hist(fMatA, "fHistA",fMeasX1,fMeasX2,fTrueX1,fTrueX2);

  // Hocker eq. (33)
  // Decompose B = QRQ' where R_{ij} = r_i^2 \delta_{ij}
  TDecompSVD QRQT(fMatB); 
  TMatrixD Q = QRQT.GetU(); // Q=U=V if B is symm. & pos.-semidef
  TVectorD rsq = QRQT.GetSig(); // R(i,i) \equiv rsq(i)

  if (0) {
    // Hocker eq. (34)
    fMatATilde = Q*fMatA;
    for (int i=0;i<fM;i++) // row index 
      for (int j=0;j<fN;j++) // col index
	fMatATilde(i,j) *= 1./TMath::Sqrt(rsq(i));
    fVecbTilde = Q*fVecb;
    for (int i=0;i<fM;i++)
      fVecbTilde(i) *= 1./TMath::Sqrt(rsq(i));
  }
  
  // Simplified scaling: assume for now cov(b) is diagonal.
  fMatATilde = fMatA;
  for (int i=0;i<fM;i++) // row index 
    for (int j=0;j<fN;j++) // col index
      fMatATilde(i,j) *= 1./TMath::Sqrt(fMatB(i,i));

  fHistATilde = Matrix2Hist(fMatATilde, "fHistATilde",
			    fMeasX1,fMeasX2,fTrueX1,fTrueX2);
  fVecbTilde = fVecb;
  for (int i=0;i<fM;i++)
    fVecbTilde(i) *= 1./TMath::Sqrt(fMatB(i,i));
  
  fHistbTilde = Vec2Hist(fVecbTilde, fMeasX1, fMeasX2, 
			 "fHistbTilde", "Scaled measured distribution");
}

void 
UnfoldingUtils::SVDAnalysis(TH2* hA, TH1* hb, TObjArray* output)
{
  // Decompose A as U*Sigma*V' and study the components. If A is m x
  // n, then U is a column-orthogonal m x n matrix, Sigma is a
  // diagonal n x n matrix (represented below as a vector), and V is
  // a column-orthonormal n x n matrix (V'*V = 1).

  TMatrixD A = Hist2Matrix(hA);
  TVectorD b = Hist2Vec(hb);

  TDecompSVD decomp(A);
  TVectorD sigma_vec = decomp.GetSig();
  TMatrixD U = decomp.GetU();
  TMatrixD VT = decomp.GetV(); VT.T();
  TMatrixD UT = U; U.T();
  int nvals = sigma_vec.GetNoElements();
  int mU = U.GetNrows(), nU = U.GetNcols(), 
    mVT = VT.GetNrows(), nVT = VT.GetNcols();
  Printf("A (%d x %d) = U (%d x %d) * Sigma (%d x %d) * V^T (%d x %d)", 
	 A.GetNrows(), A.GetNcols(), mU, nU, nvals, nvals, mVT, nVT);

  TH1D* hSigma = Vec2Hist(sigma_vec, 0., (double)nvals, "hSigma", "#sigma_{i} ");
  TH1D* hUb    = new TH1D("hUb", "|u^{T}_{i}*b| ", mU, 0, mU);
  TH1D* hUbSig = new TH1D("hUbSig", "|u^{T}_{i}*b| / #sigma_{i} ", mU, 0, mU);

  SetTH1Props(hSigma, kBlack, 0, kBlack, kFullCircle, 1.0);
  SetTH1Props(hUb, kBlue, 0, kBlue, kFullSquare, 1.0);
  SetTH1Props(hUbSig, kRed, 0, kRed, kOpenSquare, 1.0);

  output->Add(hSigma);
  output->Add(hUb);
  output->Add(hUbSig);

  for (int i=0; i<nU; i++) {
    TVectorD ui = TMatrixDColumn(U, i);
    TH1D* hui = Vec2Hist(ui, 0., (double)mU, Form("SV_u_%d",i), Form("SV_u_{%d}",i));
    output->Add(hui);
    
    // |u'_i*b|
    TVectorD uTi = TMatrixDColumn(UT, i);
    double val = TMath::Abs(uTi*b);
    double sig = hSigma->GetBinContent(i+1);
    double r = sig ? val/sig : 0.;
    hUb->SetBinContent(i+1, val);
    hUbSig->SetBinContent(i+1, r);
  }

  return;
}

TCanvas* 
UnfoldingUtils::DrawSVDPlot(TObjArray* svdhists, double ymin, double ymax)
{
  static int i=0; i++;
  TCanvas* c = new TCanvas(Form("csvd%d",i), Form("csvd%d",i), 1);
  TH1D* hSigma = (TH1D*)svdhists->At(0);//FindObject("hSigma");
  TH1D* hUb    = (TH1D*)svdhists->At(1);//FindObject("hUb");
  TH1D* hUbSig = (TH1D*)svdhists->At(2);//FindObject("hUbSig");
  if (!hSigma)
    Warning("DrawSVDPlot()","!hSigma");
  if (!hUb)
    Warning("DrawSVDPlot()","!hUb");
  if (!hUbSig)
    Warning("DrawSVDPlot()","!hUbSig");
  if (hSigma) {
    TH2F* hsvd = new TH2F(Form("hsvd%d",i), "Response matrix Picard plot;column index i;", 
			  200, 0, hSigma->GetNbinsX(), 200, ymin, ymax);
    hsvd->Draw();
    hSigma->Draw("plsame");
    hUb->Draw("plsame");
    hUbSig->Draw("plsame");
    gPad->SetLogy();
    
    TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
    leg->AddEntry(hSigma, hSigma->GetTitle(), "p");
    leg->AddEntry(hUb, hUb->GetTitle(), "ep");
    leg->AddEntry(hUbSig, hUbSig->GetTitle(), "ep");
    leg->SetFillColor(kNone);
    leg->Draw();
  }
  return c;
}

TCanvas* 
UnfoldingUtils::DrawGSVDPlot(TObjArray* svdhists, double ymin, double ymax)
{
  static int i=0; i++;
  TCanvas* c = new TCanvas(Form("cgsvd%d",i), Form("cgsvd%d",i), 1);
  TH1D* hs = (TH1D*)svdhists->At(0); // Sing. values
  TH1D* hd = (TH1D*)svdhists->At(1); // d (lambda=0)
  TH1D* hl = (TH1D*)svdhists->At(2); // d (lambda)
  if (!hs) Warning("DrawGSVDPlot()","!hs");
  if (!hd) Warning("DrawGSVDPlot()","!hd");
  if (!hl) Warning("DrawGSVDPlot()","!hl");
  if (hs) {
    TH2F* hsvd = new TH2F(Form("hgsvd%d",i), "GSVD components;column index i;", 
			  200, 0, hs->GetNbinsX(), 200, ymin, ymax);
    hsvd->Draw();
    hs->Draw("plsame");
    hd->Draw("plsame");
    hl->Draw("plsame");
    gPad->SetLogy();
    
    TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
    leg->AddEntry(hs, hs->GetTitle(), "p");
    leg->AddEntry(hd, hd->GetTitle(), "ep");
    leg->AddEntry(hl, hl->GetTitle(), "ep");
    leg->SetFillColor(kNone);
    leg->Draw();
  }
  return c;
}

TObjArray*
UnfoldingUtils::ShawSystem(const int n, double noise)
{
  TObjArray* arr = new TObjArray();
  TMatrixD A(n,n);
  TVectorD x(n);
  TVectorD b_ideal(n);
  ShawSystem(n, A, x, b_ideal, 0.);
  TVectorD b(b_ideal);

  // Add Gaussian white noise to b
  TRandom3 r3;
  for (int j=0; j<n; j++) b(j) += noise*r3.Gaus();
  
  TH2D* hA = Matrix2Hist(A, "Shaw_A",0.,1.,0.,1.);
  TH1D* hx = Vec2Hist(x, 0., 1., "Shaw_x","Truth x fn.");
  TH1D* hb = Vec2Hist(b, 0., 1., "Shaw_b","Meas. b fn.");
  TH1D* hi = Vec2Hist(b_ideal, 0., 1., "Shaw_b_ideal","Meas. b fn.");

  for (int j=0; j<n; j++)
    hb->SetBinError(j+1, noise);
  
  arr->Add(hA);
  arr->Add(hx);
  arr->Add(hb);
  arr->Add(hi);

  return arr;
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
  TRandom3 r3;
  for (int j=0; j<n; j++) b(j) += noise*r3.Gaus();

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
  
  // Unmodified chi^2 (Ax-b)'*B*(Ax-b) (Hocker eq. (30))
  // TODO: add option to select A, Ahat, Atilde, etc.
  TVectorD resid = fMatATilde*x - fVecbTilde;
  //  resid *= fMatB;
  chi2 = resid*resid;
  
  // Additive chi^2 modifier (reg. penalty value)  
  double mod = fRegWeight*SmoothingNorm(x, fRegType);
  
  chi2 += mod*mod;
  return chi2;
}

TH1D*
UnfoldingUtils::UnfoldChiSqMin(TH2* hA, TH1* hb, TH1* hXStart, TH1* hEff, TH1* hXini, 
			       double regWt, TObjArray* /*output*/, TString /*opt*/)
{
  static int id=0; id++;

  fRegWeight = regWt;
  Info("UnfoldingUtils::UnfoldChiSqMin()",
       "Using reg type %d, weight %g",fRegType, fRegWeight);

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

  // TODO: Provide option to normalize hXStart here?
  //

  // TODO: Extract min val from ICs?
  // 

  // Set up the chi^2 minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.00001);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(this, &UnfoldingUtils::RegChi2, nPars);
  min->SetFunction(f);
  
  // Initialize tmx array and pass to minimizer.
  double* tmx = new double[nBinsT+1]; // extra for overflow (?????)
  double stepSize = 0.1;
  for (int i=0; i<nBinsT; i++) {
    tmx[i] = hXStart->GetBinContent(i+1);

    if (hXini)
      tmx[i] /= hXini->GetBinContent(i+1);
    
    // Require all parameters to have a minimum positive value
    if (tmx[i] < 0) {
      tmx[i] = -tmx[i];
    }

    min->SetVariable(i, Form("xpar%d",i), tmx[i], stepSize);

    // TODO: Add case for inefficiency here?
    // 
  }

  Info("UnfoldingUtils::UnfoldChiSqMin()", 
       "Initial (regularized) chi squared = %g", RegChi2(tmx));
  min->Minimize(); 
  
  for (int i=0; i<nPars; i++) {
    tmx[i] = min->X()[i];
    double val = tmx[i];
    double err = 2*min->Errors()[i]*tmx[i];
    if (hEff) {
      double e = hEff->GetBinContent(i+1);
      val *= (e>0)? 1./e : 0.;
      err *= (e>0)? 1./e : 0.;
    }
    hUnf->SetBinContent(i+1, val);
    hUnf->SetBinError(i+1, err);
  }
  TVectorD x = Hist2Vec(hUnf);
  
  if (hXini)
    hUnf->Multiply(hXini);
  
  Info("UnfoldingUtils::UnfoldChiSqMin()", 
       "Final (regularized) chi squared = %g, curvature = %g",
       RegChi2(tmx), Curvature(x));
  return hUnf;
}

TH1D* 
UnfoldingUtils::UnfoldSVD(TH2* hA, TH1* hb, TH1* hXini, 
			  double lambda, TObjArray* output, TString opt)
{
  static int id=0; id++; // So this fn. can be called more than once
  
  int matrixType = k2DerivBCR; // favor reflected w at boundaries
  if (opt.Contains("BC0"))
    matrixType = k2DerivBC0;   // favor w=0 at boundaries
  
  // Setup
  TMatrixD A       = Hist2Matrix(hA);
  TVectorD b       = Hist2Vec(hb);
  TMatrixD L       = LMatrix(A.GetNcols(), matrixType, 1e-5);
  TMatrixD Linv    = MoorePenroseInverse(L);
  TMatrixD LTi(Linv); LTi.T(); // L^{-T}
  
  TVectorD xini(A.GetNcols());
  if (hXini) 
    xini = Hist2Vec(hXini);
  
  TVectorD eb(b);
  for (int i=0; i<b.GetNrows(); i++) 
    eb(i) = hb->GetBinError(i+i);

  //  double factor = b.Sum()/A.Sum();
  //  double factor = b.Sum()/A.Sum()/xini.Sum();
  // xini *= factor;
  //  A *= factor;
  
  if (hXini)  
    A = MultRowsByVector(A, xini);
  
  // Scale by b uncertainty (A, b --> Atilde, btilde)  
  A = DivColsByVector(A, eb);
  b = ElemDiv(b,eb);    // b_i / eb_i
  eb = ElemDiv(eb,eb);    // b_i / eb_i
  
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

  // Abs. value vectors - for analysis
  TVectorD absd(nd); // absd(i) is |(U^T*b)_i|
  TVectorD dlam(nd); // dlam(i) is tf_i*|d_i|, or s_i*|z_i|

  // Compute filter factors and z
  double s0 = s_vec(0);
  double si = 0;
  double smin = 1e-12*s0;
  for (int i=0; i<nd; i++) {
    si = s_vec(i);
    if (si < smin) si = smin;
    tf(i) = si*si/(si*si + lambda*lambda);
    z(i) = d(i)/si * tf(i);
    
    // Extra - not part of solution
    absd(i) = TMath::Abs(d(i));
    dlam(i) = absd(i)*tf(i);
  }

  // absd *= b.Sum()/xini.Sum();
  // dlam *= b.Sum()/xini.Sum();

  // Compute final solutions
  TVectorD w = Linv * V * z;
  TVectorD x(w);
  if (hXini) 
    x = ElemMult(xini,w);

  // and covariance matrices: 
  TMatrixD Z(nd,nd); 
  for (int i=0; i<nd; i++) Z(i,i) = tf(i);
  TMatrixD Wcov = Linv*V*Z*VT*LTi;
  // TVectorD wx = Wcov*xini;
  // TMatrixD Xcov = OuterProduct(xini, wx);

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
    TH1D* hd = Vec2Hist(absd,  0., nd, Form("hd%d",id), Form("|d_{i}| = |(U^{T}b)_{i}|" ));
    TH1D* hl = Vec2Hist(dlam,  0., nd, Form("hl%d",id), Form("|d^{(#lambda)}_{i}|, #lambda = %g ",lambda));
    TH1D* hw = Vec2Hist(w,     0., nd, Form("hw%d",id), Form("w^{#lambda = %g} ",lambda));
    TH1D* ht = Vec2Hist(tf,    0., nd, Form("ht%d",id), Form("s_{i}^{2}/(s_{i}^{2}+#lambda^{2}), #lambda = %g ",lambda));
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
    output->Add(hWcov);
    output->Add(hXcov);
    output->Add(hXinv);
  }
  
  TH1D* hx = Vec2Hist(x,fTrueX1,fTrueX2,Form("hSVD%d",id), Form("#lambda = %g", lambda));
  return hx;
}

TH1D* 
UnfoldingUtils::UnfoldRichardsonLucy(const TH2* hResp,
				     const TH1* hMeas,
				     const TH1* hXStart, 
				     const int nIterations, 
				     TObjArray* hists,
				     TObjArray* extras, 
				     const TH1* hXini)
{
  // See eq. (4.3), J. Bardsley & C. Vogel, 
  // SIAM J. SCI. COMPUT. Vol. 25, No. 4, pp. 1326–1343
  
  static int ihist = 0; ihist++;       // Unique ID

  TMatrixD A  = Hist2Matrix(hResp);
  TVectorD b  = Hist2Vec(hMeas);
  TVectorD x0 = Hist2Vec(hXStart);
  double x1 = hXStart->GetXaxis()->GetXmin();
  double x2 = hXStart->GetXaxis()->GetXmax();
  int m = A.GetNrows(); // # meas bins
  int n = A.GetNcols(); // # true bins
  TMatrixD AT(A); AT.T();
  TVectorD x(x0); // Must be positive
  TVectorD xini(n);
  if (hXini) {
    xini = Hist2Vec(hXini);
    A = MultRowsByVector(A, xini);
  }
  TVectorD eb(b);
  for (int i=0; i<b.GetNrows(); i++) 
    eb(i) = hMeas->GetBinError(i+i);

  // If requested, divide by b uncertainty (A, b --> Atilde, btilde)  
  //  if (opt.Contains("SB")) {
  if (0) { // This does not work
    A = DivColsByVector(A, eb);
    b = ElemDiv(b,eb);
    eb = ElemDiv(eb,eb);
  }
  TVectorD bvar = ElemMult(eb,eb);  
  b += bvar;
  TVectorD ones(m); 
  for (int i=0; i<m; i++) 
    ones(i)=1.0;

  // Construct d vector whose elements are the column sums in A
  TVectorD d(n);   // sum over rows
  for (int i=0; i<m; i++) {   // row loop
    for (int j=0; j<n; j++) { // col loop
      d(j) += A(i,j);
    }
  }

  // L-curve graph - one point per iteration
  TGraph* gL = 0;
  if (extras) {
    gL = new TGraph(nIterations);
    gL->SetLineColor(kRed);
    gL->SetMarkerColor(kRed);
    gL->SetMarkerStyle(kFullCircle);
    gL->SetMarkerSize(1.5);
    gL->SetLineWidth(2);
    gL->SetNameTitle(Form("LCurve_RL_%d",ihist), 
		     Form("Richardson-Lucy L-Curve;"
			  "Residual norm ||Ax_{k}-b||_{2};"
			  "Solution norm ||x_{k}||_{2}"));
  }

  for (int k=0; k<=nIterations; k++) {
    printf("Richardson-Lucy iteration %d\r", k);
  
    if (k > 0) {
      TVectorD Ax = A*x + bvar; 
      TVectorD r = ElemDiv(b, Ax);
      TVectorD ATr = AT*r;
      TVectorD AT1 = AT*ones;
      TVectorD xoverAT1 = ElemDiv(x, AT1);
      x = ElemMult(xoverAT1, ATr);

      if (gL) {    // Set L-Curve points
	TVectorD xrescaled = ElemMult(x,xini); // only nonzero if hXini exists
	TVectorD resid = (hXini) ? A*xrescaled-b : A*x-b;
	gL->SetPoint(k-1, TMath::Sqrt(resid*resid), TMath::Sqrt(x*x));
      }
    }

    if (hists) { // Save nth iterate as a TH1D
      if (hXini) {
	TVectorD xrescaled = ElemMult(x,xini);
	hists->Add(XHist(xrescaled,Form("RL%d",ihist),k,x1,x2,0.,""));
      }
      else
	hists->Add(XHist(x,Form("RL%d",ihist),k,x1,x2,0.,""));
    }
  } // end iteration loop
  cout << endl;

  TH1D* hx=0;
  if (hists)
    hx = (TH1D*)hists->At(nIterations);
  else {
    if (hXini) {
      TVectorD xrescaled = ElemMult(x,xini);
      hx = XHist(xrescaled,Form("RL%d",ihist),nIterations,x1,x2,0.,"");
    }
    else
      hx = XHist(x,Form("RL%d",ihist),nIterations,x1,x2,0.,"");
  }
  if (extras) {
    if (gL) extras->Add(gL);
  }

  return hx;
}

TH1D* 
UnfoldingUtils::UnfoldPCGLS(const TH2* hA, 
			    const TH1* hb, 
			    TH1* hXini,
			    const int nIterations, 
			    TObjArray* hists,
			    TObjArray* extras, 
			    int LType,
			    TString opt)
{
  //
  // Key ingredients and their dimensions
  // -----------------------------------------------------------------
  // A: m x n  (hResp) TH2 is true (y) vs meas (x), as in RooUnfold.
  // b: m      (hMeas)
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
  TMatrixD A       = Hist2Matrix(hA);
  TVectorD b       = Hist2Vec(hb);
  TMatrixD L       = LMatrix(A.GetNcols(), LType, 1e-5);
  TMatrixD W       = Null(L);          // Nullspace of L
  TVectorD eb(b);
  for (int i=0; i<b.GetNrows(); i++) 
    eb(i) = hb->GetBinError(i+i);
  
  int ncols = A.GetNcols();
  int p = L.GetNrows();
  double x1 = fTrueX1, x2 = fTrueX2;
  TVectorD xini(ncols);
  if (hXini) {
    x1 = hXini->GetXaxis()->GetXmin();
    x2 = hXini->GetXaxis()->GetXmax();
    xini = Hist2Vec(hXini);
    A = MultRowsByVector(A, xini);
  }

  // If requested, divide by b uncertainty (A, b --> Atilde, btilde)  
  if (opt.Contains("SB")) {
    A = DivColsByVector(A, eb);
    b = ElemDiv(b,eb);
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
			  "Solution norm ||x_{k}||_{2}"));
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
	gL->SetPoint(n-1, TMath::Sqrt(r*r), TMath::Sqrt(x*x));
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
    TMatrixD Linv(L);
    Linv.Invert();
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
UnfoldingUtils::ElemDiv(const TVectorD& x, const TVectorD& y)
{
  int nx = x.GetNoElements();
  int ny = y.GetNoElements();
  if (nx != ny) {
    gROOT->Error("ElemDiv()", "mismatch nx=%d, ny=%d", nx, ny);
    gSystem->Exit(-1);
  }
  TVectorD result(nx);
  for (int i=0; i<nx; i++) {
    result(i) = (y(i) > 1e-15) ? x(i) / y(i) : 0.;
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
  // The factor normto could be hMeas->Integral();
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

TMatrixD 
UnfoldingUtils::Null(TMatrixD& A)
{
  // Return matrix whose columns form the orthonormal basis for the
  // nullspace of A. Perform SVD on A so that A = USV'. The columns of
  // V whose corresponding singular values are "zero" are the null basis
  // vectors.
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
UnfoldingUtils::BandedDiagonalMatrix(TH1* hDpt, const int nMeas, const int nTrue)
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
  TH2D* h = Matrix2Hist(Abkg, "banded_diag", 0, wdpt*nMeas, 0, wdpt*nTrue); 
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
