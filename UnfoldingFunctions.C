// NOTE: This code is deprecated. It has been replaced by the UnfoldingUtils source files.

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TNamed.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TDecompLU.h"
#include "TDecompQRH.h"
#include "TGraphTime.h"
#include "TSystem.h"

#include <iostream>
using namespace std;

// Function prototypes
TMatrixD Hist2Matrix(const TH2* h);
TMatrixD Hist2Matrix(const TH2* h, Int_t nx, Int_t ny, const TH1* norm, Bool_t overflow = kFALSE);
void Vec2Hist(const TVectorD& v, TH1* h, Int_t nb, Bool_t overflow = kFALSE);
TVectorD Hist2Vec(const TH1* h, Int_t nb, Bool_t overflow = kFALSE);
double GetBinContent(const TH1* h, Int_t i, Bool_t overflow = kFALSE);
int GetBin(const TH1* h, Int_t i, Bool_t overflow = kFALSE);
int GetBinDim(const TH1* h, Int_t i);

// TGraphTime* Animation(TObjArray* moveObjs, TObjArray* statObjs, 
// 		      TString opt="");
TMatrixD MoorePenroseInverse(TMatrixD& A, double tol = 1e-15);
TMatrixD* GetPseudoInverse(TMatrixD& A);
TMatrixD* GetDerivativeMatrix(int n, int d, TString opt="");
TMatrixD* GetNullSpaceMatrix(int n, int d);
void LTSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y);
void LSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y, 
	    const TMatrixD& W, const TMatrixD& T);
TVectorD ElemMult(const TVectorD& x, const TVectorD& y);
TVectorD ElemDiv (const TVectorD& x, const TVectorD& y);

double UnfoldingScaleFactor(TH2* h, int definition = 0);
TH2D* Transpose(TH2* h);
TH1D* InnerProduct(TH2* hR, TH1* hT, TH1* hResult = 0);
TH2D* MatrixProduct(TH2* hA, TH2* hB, TH2* hResult = 0);
void SetStatError(TH1* h, int nEntries);
TH1D* TrueInputHist(int nbins, double xmin, double xmax, 
		    TFormula form, int nEnt, TString opt="stat");
TH2D* ResponseMatrix(int nx, double xmin, double xmax, 
		     int ny, double ymin, double ymax, 
		     double norm, double sigma, double ptfac, 
		     TString opt);
TH1* UnfoldPCGLS(const TMatrixD& A, const TVectorD& b, 
		 const int nIterations, TObjArray* hists, 
		 TObjArray* extras, const TMatrixD& L, 
		 const TMatrixD& W, const double x1, 
		 const double x2, const double normTo = 0.);
TH1* UnfoldRichardsonLucy(const TMatrixD& A, 
			  const TVectorD& b, 
			  const TVectorD& x0, 
			  const int nIterations, 
			  TObjArray* hists,
			  TObjArray* extras, 
			  const double x1, 
			  const double x2,
			  const double normTo = 0.);
// TH1* UnfoldTSVD(const TMatrixD& A, const TVectorD& b, 
// 		const int nIterations, TObjArray* hists, 
// 		TObjArray* extras);
TGraph* LCurve(TMatrixD& A, TVectorD& b, TObjArray* xhists, TString name);
TGraph* L2NormVsStep(TH1* hTrue, TObjArray* xhists, TString opt="");
TH1D* XHist(TVectorD& x, TString base, int k, double xMin, 
	    double xMax, double normto, TString opt="");
void     HistFromVector(const TVectorD& v, TH1* h);
void     ArrayFromTH1(const TH1* h, double val[], double err[]);
void     ArrayFromTH2(const TH2* h, double val[], double err[]);
void     MatrixFromTH2(const TH2& h, TMatrixD& A);
void     VectorFromTH1(const TH1& h, TVectorD& v);

void SetHistProps(TH1* h,
		  Int_t linecolor,
		  Int_t fillcolor,
		  Int_t markercolor,
		  Int_t markerstyle,
		  Double_t markersize);



void GetRangeTH2(const TH2* h, 
		 int& nx, double& x1, double& x2, 
		 int& ny, double& y1, double& y2)
{
  nx = h->GetNbinsX();
  ny = h->GetNbinsY();
  x1 = h->GetXaxis()->GetXmin();
  x2 = h->GetXaxis()->GetXmax();
  y1 = h->GetYaxis()->GetXmin();
  y2 = h->GetYaxis()->GetXmax();
}

TGraph* LCurve(TMatrixD& A, TVectorD& b, TObjArray* xhists, TString name)
{
  TGraph* g = new TGraph();
  static int ihist = 0; ihist++;       // Unique ID
  for (int i=0; i<xhists->GetEntries(); i++) {
    TH1* h = (TH1*) xhists->At(i);
    int n = h->GetNbinsX();
    //    TVectorD* x = RooUnfoldResponse::H2V(h, n, kFALSE);
    TVectorD x = Hist2Vec(h, n, kFALSE);
    
    int m = A.GetNrows();
    int nb = b.GetNoElements();
    if (m != nb)
      gROOT->Error("LCurve()", "%d rows in A != %d", m, nb);

    TVectorD r = A*x-b;
    g->SetPoint(i, TMath::Sqrt(r*r), TMath::Sqrt(x*x));
  }
  g->SetLineColor(kBlue);
  g->SetMarkerColor(kBlue);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(2);
  g->SetNameTitle(Form("LCurve_%s_%d",name.Data(), ihist), 
		  Form("%s L-Curve;"
		       "Residual norm ||Ax_{k}-b||_{2};"
		       "Solution norm ||x_{k}||_{2}", name.Data()));
  return g;
}

TVectorD ElemDiv(const TVectorD& x, const TVectorD& y)
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

TVectorD ElemMult(const TVectorD& x, const TVectorD& y)
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

// TObjArray* ComputeRescaledSystem(TH2* hAhat, TH1* hXini, TH1* hMeas)
// {
//   // hAhat must follow RooUnfold convention for TH2: x=meas, y=true
//   int m = hAhat->GetNbinsX();     // Measured bins (rows)
//   int n = hAhat->GetNbinsY();     // True/gen bins (cols)
//   int nXini = hXini->GetNbinsX(); // Must equal n
//   int nMeas = hMeas->GetNbinsX(); // Must equal m

//   // Compute B = cov(b) assuming only indep. stat. err. at this time
//   TMatrixD B(m,m); B.UnitMatrix();
//   for (int i=0;i<m;i++) 
//     B(i,i) *= hMeas->GetBinError(i+1)*hMeas->GetBinError(i+1);

//   // Scale the probability response matrix Ahat so columns sum to counts
//   // (instead of 1.0). This "cancels" the scaling w_j = x_j/xini_j;
//   // see A. Hocker, NIM A 372 (1996) 469-481.
//   TMatrixD Ahat = Hist2Matrix(hAhat);
//   TMatrixD A(Ahat);
//   for (int i=0;i<m;i++) // row index 
//     for (int j=0;j<n;j++) // col index
//       A(i,j) *= hXini->GetBinContent(j+1);

//   // Hocker eq. (33)
//   // Decompose B = QRQ' where R_{ij} = r_i^2 \delta_{ij}
//   TDecompSVD QRQT(B); 
//   TMatrixD Q = QRQT.GetU(); // Q=U=V if B is symm. & pos.-semidef
//   TVectorD rsq = QRQT.GetSig(); // R(i,i) \equiv rsq(i)

//   // Hocker eq. (34)
//   TMatrixD Atilde(Q*Ai);
//   for (int i=0;i<m;i++) // row index 
//     for (int j=0;j<n;j++) // col index
//       Atilde(i,j) *= 1./TMath::Sqrt(rsq(i));
//   TVectorD btilde = Q*b;
//   for (int i=0;i<m;i++)
//     btilde(i) *= 1./TMath::Sqrt(rsq(i));

//   TObjArray* results = new TObjArray();
//   results->Add(&Ahat);    // 0 (m)
//   results->Add(&A);       // 1 (m)
//   results->Add(&Atilde);  // 2 (m)
//   results->Add(&b);       // 3 (v)
//   results->Add(&btilde);  // 4 (v)
//   results->Add(&B);       // 5 (m)

//   return results;
// }

// TMatrixD GetAhat(TObjArray* arr)
// {
//   return (TMatrixD)arr->At(0);
// }
// TMatrixD GetA(TObjArray* arr)
// {
//   return (TMatrixD)arr->At(1);
// }
// TMatrixD GetATilde(TObjArray* arr)
// {
//   return (TMatrixD)arr->At(2);
// }
// TVectorD Getb(TObjArray* arr)
// {
//   return (TVectorD)arr->At(3);
// }
// TVectorD GetbTilde(TObjArray* arr)
// {
//   return (TVectorD)arr->At(4);
// }
// TMatrixD GetB(TObjArray* arr)
// {
//   return (TMatrixD)arr->At(5);
// }

TH1* UnfoldRichardsonLucy(const TMatrixD& A, 
			  const TVectorD& b,
			  const TVectorD& x0, 
			  const int nIterations, 
			  TObjArray* hists,
			  TObjArray* extras, 
			  const double x1, 
			  const double x2,
			  const double normTo)
{
  // For k=0,1,2,..., compute x[k] as
  // x[k+1] = x[k] edot A'(b ediv (Ax[k])) ediv d,
  // where edot (ediv) is element-wise multiplication (division). 

  static int ihist = 0; ihist++;       // Unique ID
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
  int m = A.GetNrows(); // # meas bins
  int n = A.GetNcols(); // # true bins
  TMatrixD AT(A); AT.T();
  TVectorD x(x0); // Must be positive

  // // Initialize x. x(j) must be > 0.
  // for (int j=0; j<n; j++) x(j) = x0(j);

  // Construct d vector whose elements are the column sums in A
  TVectorD d(n);   // sum over rows
  for (int i=0; i<m; i++) {   // row loop
    for (int j=0; j<n; j++) { // col loop
      d(j) += A(i,j);
    }
  }

  for (int k=0; k<=nIterations; k++) {
    printf("Richardson-Lucy iteration %d\r", k);
  
    if (k > 0) {
      TVectorD Ax = A*x;
      TVectorD r = ElemDiv(b, Ax);
      TVectorD ATr = AT*r;
      TVectorD ATrd = ElemDiv(ATr, d);
      TVectorD xnew = ElemMult(x, ATrd);
      x = xnew;

      if (gL) {    // Set L-Curve points
	TVectorD resid = A*x-b; // for the L-Curve
	gL->SetPoint(k-1, TMath::Sqrt(resid*resid), TMath::Sqrt(x*x));
      }
    }

    if (hists) // Save nth iterate as a TH1D
      hists->Add(XHist(x,Form("RL%d",ihist),k,x1,x2,normTo,""));
  } // end iteration loop
  cout << endl;

  TH1D* hx=0;
  if (hists)
    hx = (TH1D*)hists->At(nIterations);
  else
    hx = XHist(x,Form("RL%d",ihist),nIterations,x1,x2,normTo,"");
  if (extras) {
    if (gL) extras->Add(gL);
  }
  return hx;
}

TH1* UnfoldPCGLS(TH1* hMeas, 
		 TH2* hResp, 
		 int nIterations, 
		 TObjArray* hists,
		 TObjArray* extras,
		 TMatrixD* Lptr,
		 TMatrixD* Wptr)
{
  int nbinsx, nbinsy;
  double x1, x2, y1, y2;
  GetRangeTH2(hResp,nbinsx,x1,x2,nbinsy,y1,y2);
  int nrows = nbinsx, ncols = nbinsy; // nrows (m) x ncols (n) of A.
  TH1* hNorm=0;                       // TODO: implement efficiency

  // Response matrix and measured distribution
  // TMatrixD* Aptr = RooUnfoldResponse::H2M(hResp, nrows, ncols, hNorm, kFALSE);
  // TVectorD* bptr = RooUnfoldResponse::H2V(hMeas, nrows, kFALSE);
  // TMatrixD A(*Aptr);
  // TVectorD b(*bptr);
  // delete Aptr; 
  // delete bptr;

  TMatrixD A = Hist2Matrix(hResp, nrows, ncols, hNorm, kFALSE);
  TVectorD b = Hist2Vec(hMeas, nrows, kFALSE);

  TMatrixD L, W;
  if (Lptr) {
    L.ResizeTo(*Lptr); L = *Lptr;
    if (Wptr) {
      W.ResizeTo(*Wptr); W = *Wptr;
    }
  }
  else {
    L.ResizeTo(ncols, ncols); 
    L.UnitMatrix();
    W.ResizeTo(0,0);
  }
  
  double normTo = hMeas->Integral() * double(ncols)/nrows;

  return UnfoldPCGLS(A, b, nIterations, hists, extras, L, W, x1, x2, normTo);
}

TH1* UnfoldPCGLS(const TMatrixD& A, 
		 const TVectorD& b, 
		 const int nIterations, 
		 TObjArray* hists,
		 TObjArray* extras, 
		 const TMatrixD& L, 
		 const TMatrixD& W, 
		 const double x1, 
		 const double x2,
		 const double normTo)
{
  //
  // Key ingredients and their dimensions
  // -----------------------------------------------------------------
  // A: m x n  (hResp) TH2 is true (y) vs meas (x), as for RooUnfold.
  // b: m      (hMeas)
  // x: n      (solutions)
  // L: p x n  p has no direct restrictions
  // W: n x ?  Basis vectors of L nullspace. W != 0 for p < n.
  // -----------------------------------------------------------------
  //
  // Some equations below reference the book "Rank Deficient and
  // Discrete Ill-Posed Problems" by P.C. Hansen, sections 2.3.2 and
  // 6.3.
  //

  bool reorth = 1;
  static int ihist = 0; ihist++;       // Unique ID
  int ncols = A.GetNcols();
  int p = L.GetNrows();
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
    T = AWinv*A;                  // eq. 2.47 (NAA in pcgls.m)
    x0 = W*AWinv*b;               // eq. 2.46
  }
  
  TVectorD x(x0);                    // kth solution
  TMatrixD AT(TMatrixD::kTransposed, A);
  TVectorD r = b - A*x0;
  TVectorD s = AT*r;
  TVectorD q(p), q1(p);

  //  cout<<"LTSolve"<<endl;  
  LTSolve(q1, L, s);        // q1 = L11^{-T} * s. length p.

  //  cout<<"LSolve"<<endl;  
  LSolve(q, L, q1, W, T);   // q = L_A^+ * q1

  TVectorD z(q);
  double dq = s*q, dq2=0;
  double alpha=0, beta=0;

  if (reorth) {
    double q1norm = TMath::Sqrt(q1*q1);
    if (q1norm > 0) q1norm = 1./q1norm;
    TMatrixDColumn(Q1n, 0) = q1norm*q1;
  }

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
      if (reorth) {
	for (int i=0; i<n; i++) {
	  TVectorD qi = TMatrixDColumn(Q1n, i);
	  q1 -= (qi*q1)*qi; // Modified Graham-Schmidt
	  double q1norm = TMath::Sqrt(q1*q1);
	  if (q1norm > 0) q1norm = 1./q1norm;
	  TMatrixDColumn(Q1n, n) = q1norm*q1;
	}
      }
      LSolve(q, L, q1, W, T);
      dq2 = s*q;
      beta = dq2/dq;
      dq = dq2;
      z = q + beta*z;

      if (gL)    // Set L-Curve points
	gL->SetPoint(n-1, TMath::Sqrt(r*r), TMath::Sqrt(x*x));
    }

    if (hists) // Save nth iterate as a TH1D
      hists->Add(XHist(x,Form("PCGLS%d",ihist),n,x1,x2,normTo,""));
  } // end iteration loop
  cout << endl;

  TH1D* hx=0;
  if (hists)
    hx = (TH1D*)hists->At(nIterations);
  else
    hx = XHist(x,Form("PCGLS%d",ihist),nIterations,x1,x2,normTo,"");

  if (extras) {
    if (gL) extras->Add(gL);
  }

  return hx;
}


TH1* UnfoldCGLS(TH1* hMeas, TH2* hResp, int nIterations, TObjArray* hists, TObjArray* extras)
{
  // hResp should true (y) vs meas (x), following the RooUnfold
  // convention.
  // Result is scaled to hMeas->Integral().

  // Unique ID in case of multiple calls
  static int ihist = 0; ihist++;
  int nbinsx, nbinsy;
  double x1, x2, y1, y2;
  GetRangeTH2(hResp,nbinsx,x1,x2,nbinsy,y1,y2);
  int nrows = nbinsx, ncols = nbinsy;
  TH1* hNorm=0; // TODO: implement passing efficiency into this fn.
  double normTo = hMeas->Integral() * double(ncols)/nrows;
  TGraph* gL = 0;
  if (extras) {
    gL = new TGraph(nIterations);
    gL->SetNameTitle(Form("LCurve_CGLS_%d",ihist), Form("CGLS L-Curve;||Ax_{k}-b||_{2};||x_{k}||_{2}"));
  }
  // Response matrix and measured distributions
  TMatrixD A = Hist2Matrix(hResp, nrows, ncols, hNorm, kFALSE);
  TVectorD b = Hist2Vec(hMeas, nrows, kFALSE);

  // TMatrixD* Aptr = RooUnfoldResponse::H2M(hResp, nrows, ncols, hNorm, kFALSE);
  // TVectorD* bptr = RooUnfoldResponse::H2V(hMeas, nrows, kFALSE);
  // TMatrixD A(*Aptr);
  // TVectorD b(*bptr);
  // delete Aptr;
  // delete bptr;

  TMatrixD AT(TMatrixD::kTransposed, A);
  TVectorD x(ncols);  x.Zero();       // x0 = null to start
  TVectorD r = b - A*x;
  TVectorD d = AT*r;
  TVectorD s(d);
  double normr2 = d*d;
  double normr2_new = 0;
  double alpha=0, beta=0;             // Ratios of 2-norms

  bool reorth = 1;
  TMatrixD ATr(ncols, nIterations+1); // Store d in ATr columns
  if (reorth) {
    double dnorm = TMath::Sqrt(d*d);
    if (dnorm > 0) dnorm = 1./dnorm;
    TMatrixDColumn(ATr, 0) = dnorm*d;
    gROOT->Info("CGLS()", "MGS Reorthogonalization enabled");
  }

  for (int n=0; n<=nIterations; n++) {

    if (n>0) { // Run CGLS algorithm for n=1,2,...

      // Update x and r vectors
      TVectorD Ad = A*d;
      alpha = normr2 / (Ad*Ad);
      x += alpha*d;
      r -= alpha*Ad;
      s = AT*r;

      if (reorth) { // Reorthogonalize s to previous s vectors
     	for (int i=0; i<n; i++) {
     	  TVectorD si = TMatrixDColumn(ATr, i);
     	  s -= (si*s)*si;
     	  double snorm = TMath::Sqrt(s*s);
     	  if (snorm > 0) snorm = 1./snorm;
     	  TMatrixDColumn(ATr, n) = snorm*s;
     	}
      }
      
      // Update d vector
      normr2_new = s*s;
      beta = normr2_new / normr2;
      normr2 = normr2_new;
      d = s + beta*d;

      if (gL) {
	gL->SetPoint(n-1, TMath::Sqrt(r*r), TMath::Sqrt(x*x));
      } 
    }
    if (hists) {
      hists->Add(XHist(x,Form("CGLS%d",ihist),n,x1,x2,normTo,""));
    }
  } // end iteration loop
  
  TH1D* hx=0;
  if (hists)
    hx = (TH1D*)hists->At(nIterations);
  else
    hx = XHist(x,Form("CGLS%d",ihist),nIterations,x1,x2,normTo,"");
  
  if (extras) {
    if (gL) extras->Add(gL);
  }
  
  return hx;
}

TH1D* XHist(TVectorD& x, TString base, int k, double xMin, double xMax,
	    double normto, TString /*opt*/)
{
  // base should have a form like "CGLS"+"iHist"
  // k is the iteration step
  // The factor normto could be hMeas->Integral();
  const char* name  = Form("h%s_%d", base.Data(), k);
  const char* title = Form("%s (%d iterations)", base.Data(), k);
  TH1D* hx = new TH1D(name, title, x.GetNoElements(), xMin, xMax);
  HistFromVector(x, hx);
  if (normto) {
    double norm = hx->Integral() ? normto/hx->Integral() : 0;
    hx->Scale(norm);
  }
  // TODO assign proper stat uncertainty
  return hx;
}

TH1* UnfoldLandweber(const TH1* hMeas, const TH2* hResp, const TH1* hNorm,
		     const double omega, const int nIterations, 
		     TObjArray* hists,  int sampleInterval,
		     TString opt)
{
  // hResp should true (y) vs meas (x), following the RooUnfold
  // convention.

  static int ihist = 1;
  int nbinsx, nbinsy;
  double x1, x2, y1, y2;
  GetRangeTH2(hResp,nbinsx,x1,x2,nbinsy,y1,y2);
  int nrows = nbinsx, ncols = nbinsy;
  double normTo = hMeas->Integral() * double(ncols)/nrows;

  if (nrows != hMeas->GetNbinsX()) { // b must have nbinsy
    gROOT->Error("UnfoldLandweber()",
		 "%d y bins in %s != %d in %s",
		 nbinsy, hResp->GetName(), hMeas->GetNbinsX(), hMeas->GetName());
    return 0;
  }

  // Response matrix and measured distributions
  TMatrixD A = Hist2Matrix(hResp, nrows, ncols, hNorm, kFALSE);
  TVectorD b = Hist2Vec(hMeas, nrows, kFALSE);
  // TMatrixD* A = RooUnfoldResponse::H2M(hResp, nrows, ncols, hNorm, kFALSE);
  // TVectorD* b = RooUnfoldResponse::H2V(hMeas, nrows, kFALSE);

  TMatrixD AT(TMatrixD::kTransposed, A);
  TMatrixD M = AT*A;

  // zeroth iterate x0 = omega * AT * b
  TVectorD x0 = omega*AT*b;
  TVectorD xn(x0);                            // solution
  TMatrixD I(TMatrixD::kUnit, M);             // identity matrix

  // Precompute M as 1 - omega*AT*A
  M = I - omega*M;

  // Compute the iteration
  // xn = xn + omega*AT*(b - A*xn)
  for (int n=0; n<=nIterations; n++) {
    
    if (n>0) {
      xn = x0 + M*xn;

      // Enforce non-negativy? E.g. for PDFs
      if (opt.Contains("+")) {
	for (int j=0; j<xn.GetNoElements(); j++)
	  if (xn(j) < 0)
	    xn(j) = 0;
      }
    }

    if (hists && (n%sampleInterval==0))
      hists->Add(XHist(xn,Form("Landweber%d",ihist),n,x1,x2,normTo,""));
    
  }

  TH1D* hx=0;
  if (hists)
    hx = (TH1D*)hists->At(nIterations);
  else
    hx = XHist(xn,Form("Landweber%d",ihist),nIterations,x1,x2,normTo,"");

  return hx;
}

TH1D* TrueInputHist(int nbins, double xmin, double xmax, TFormula formula, int nEntries, TString opt)
{
  // Generate input truth spectrum.
  
  static int ihist = 1;
  TString name = Form("hTrue%d", ihist); ihist++;
  TH1D* h = new TH1D(name.Data(), "True input spectrum", nbins, xmin, xmax);

  //  double n = 6.6, T = 0.145;
  //  double n = 6, T = 0.3;

  for (int i=1; i<=nbins; ++i) {
    double x = h->GetBinCenter(i);
    double y = formula.Eval(x);
    // if (x > 0.)
    //   y = TMath::Abs(x)*(n-1.0)*(n-2.0)/((n*T)*(n*T)) * 
    // 	TMath::Power(1.0+TMath::Abs(x)/(n*T), -n);
    
    h->SetBinContent(i, y);
  }

  if (opt.Contains("stat"))
    SetStatError(h, nEntries);
  
  return h;
}

TH1* HistWithNoise(TH1* hIn, TString /*opt*/)
{
  // Add Gaussian fluctuations to each bin of hIn, making the entries
  // integer counts.
  TH1* h = (TH1*)hIn->Clone(Form("%sWithPoissonNoise", hIn->GetName()));
  int nbins = h->GetNbinsX();
  TRandom3 rndm;
  TDatime t;
  rndm.SetSeed(t.GetTime());

  for (int i=1; i<=nbins; ++i ) {
    double val = h->GetBinContent(i);
    double err = h->GetBinError(i);
    val += rndm.Gaus(0, err);
    h->SetBinContent(i, TMath::Nint(val));
    double err2 = h->GetBinContent(i) ? TMath::Sqrt(h->GetBinContent(i)) : 1.;
    //    h->SetBinError(i, TMath::Sqrt(h->GetBinContent(i)));
    h->SetBinError(i, err2);
  }
  return h;
}

void SetStatError(TH1* h, int nEntries)
{
  int nbins = h->GetNbinsX();
  if (nEntries)
    h->Scale(double(nEntries)/h->Integral());
  for (int i=1; i<=nbins; ++i ) {
    h->SetBinError(i, TMath::Sqrt(h->GetBinContent(i)));
  }
}

// TODO: not yet working; fix
TH1D* FoldedByFunction(TH1* hTrue, TFormula fn)
{
  static int ihist = 1;
  TString name = Form("hFoldedByFunction%d", ihist); ihist++;
  int nbx = hTrue->GetNbinsX();
  int nbinsout = nbx; // could generalize later

  TH1D* h = dynamic_cast<TH1D*>(hTrue->Clone(name.Data()));
  h->Reset();

  for (int j=1; j<=nbinsout; j++) {
    double sum   = 0.0;  // convolution sum for bin j
    for (int i=1; i<=nbx; i++) {
      double val = hTrue->GetBinContent(i);
      double xi  = hTrue->GetBinLowEdge(i);
      double dx  = hTrue->GetBinWidth(i);

      //      printf("%3.1g, %.0f, %3.1g   ", val, xi, fn.Eval(xi));
      sum += val*(fn.Eval(xi) - fn.Eval(xi-dx));
  }
    h->SetBinContent(j, sum);
    cout << sum << endl;
  }
  return h;
}

TH1D* FoldedByGaussian(TH1* hTrue, double sigma)
{
  // Function modified from a C version written by
  // A. Laszlo in libunfold:
  // http://www.rmki.kfki.hu/~laszloa/
  static int ihist = 1;
  TString name = Form("hFolded%d", ihist); ihist++;
  int nbinsin = hTrue->GetNbinsX();
  int nbinsout = nbinsin;
  double xi=0.0, xj=0.0, xjmin=0.0, xjmax=0.0, Di=0.0, Dj=0.0;
  TH1D* h = dynamic_cast<TH1D*>(hTrue->Clone(name.Data()));
  h->Reset();

  // Loop over h bins (to fill them)
  for (int j=1; j<=nbinsout; ++j) {
    Dj    = h->GetBinWidth(j);
    xj    = h->GetBinCenter(j);
    xjmin = h->GetBinLowEdge(j);
    xjmax = xjmin + Dj;
    double sum   = 0.0;  // convolution sum for bin j

    // Loop over hTrue bins
    for (int i=1; i<=nbinsin; ++i) {
      Di = hTrue->GetBinWidth(i);
      xi = hTrue->GetBinCenter(i); // Gaussian mean.
      double val = hTrue->GetBinContent(i);

      // This would be the valid code if underflow/overflow is
      // counted into the histogram. It does not satisfy
      // necessary regularity conditions, so then linear
      // iterative unfolding wouldn't work.
      // foldedpdf[j]+=
      //   0.5*((j==0 ? 1.0 : erf((x-ybinmin)/(M_SQRT2*sigma)))-
      // 	   (j+1==nbinsout ? -1.0 : erf((x-ybinmax)/(M_SQRT2*sigma))))/Dy*inputpdf[i]*Dx;
	    
      /* This is the valid code if underflow/overflow is discarded. */
      sum += val*Di*0.5*(TMath::Erf((xi-xjmin)/(TMath::Sqrt2()*sigma))-
			 TMath::Erf((xi-xjmax)/(TMath::Sqrt2()*sigma)))/Dj;
    }
    h->SetBinContent(j, sum);
  }

  return h;
}

TH2D* ResponseMatrix(int nx, double xmin, double xmax, 
		     int ny, double ymin, double ymax, 
		     double norm, double s, double f, TString /*opt*/)
{
  // Create a toy response matrix ptmeas (y) vs. pttrue (x) where the
  // y axis is smeared by sqrt(s^2 + f^2 xtrue).
  static int ihist = 1;
  TString name = Form("hResp%d", ihist); ihist++;
  TH2D* h = new TH2D(name.Data(), "Response matrix",  		    
		     nx, xmin, xmax, ny, ymin, ymax);
  for (int i=1; i<=nx; ++i) { // true pt bins
    double x = h->GetXaxis()->GetBinCenter(i); // mu
    //    if (x<0) continue;
    for (int j=1; j<=ny; ++j) { // smeared pt bins
      double y = h->GetYaxis()->GetBinCenter(j);
      double z = TMath::Gaus(y, x, TMath::Sqrt(s*s + f*f*x));
      h->SetBinContent(i, j, z);
    }
  }

  h->Scale(norm/h->Integral());
  return h;
}

TH2D* ResponseMatrix(TFormula true_fn, TFormula conv_fn,
		     int nx, double xmin, double xmax, 
		     int ny, double ymin, double ymax, 
		     double norm, TString opt)
{

  // ==================== Create Response matrix =====================
  //
  // Truth axis is weighted by true_fn, measured axis is truth value
  // folded by conv_fn.
  //
  // Default: x=true, y=smeared by yfn. Choose "y=true" (a la
  // RooUnfold) for opposite convention, i.e. true_fn weights the
  // y-direction, conv_fn smears in the x-direction.
  //
  // =================================================================

  static int ihist = 1;
  TString name = Form("hResp%d", ihist); ihist++;
  TH2D* h = new TH2D(name.Data(), "Response matrix",  		    
		     nx, xmin, xmax, ny, ymin, ymax);

  double epsilon = 1e-50;

  if (opt.Contains("y=true")) {
    for (int j=1; j<=ny; ++j) {
      double y = h->GetYaxis()->GetBinCenter(j);
      for (int i=1; i<=nx; ++i) {
	double x = h->GetXaxis()->GetBinCenter(i); // mu
	double z = conv_fn.Eval(x-y);
	if (z<epsilon) z = 0;
	z *= true_fn.Eval(y); // weighting in y
	h->SetBinContent(i, j, z);
      }
    }
  }
  
  else { // default: x=true
    for (int i=1; i<=nx; ++i) { // true pt bins
      double x = h->GetXaxis()->GetBinCenter(i); // mu
      for (int j=1; j<=ny; ++j) { // smeared pt bins
	double y = h->GetYaxis()->GetBinCenter(j);
	double z = conv_fn.Eval(y-x);
	if (z<epsilon) z = 0;
	z *= true_fn.Eval(x); // weighting in x
	h->SetBinContent(i, j, z);
      }
    }
  }
  h->Scale(norm/h->Integral());
  return h;
}

/*
TH2D* Transpose(TH2*h)
{
  static int ihist = 1;
  TString name = Form("%sTranspose%d", h->GetName(), ihist); ihist++;

  TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis(); 
  double xmin = ax->GetXmin(), xmax = ax->GetXmax(); 
  double ymin = ay->GetXmin(), ymax = ay->GetXmax(); 
  int nx = h->GetNbinsX(), ny = h->GetNbinsY();
  TH2D* ht = new TH2D(name.Data(), name.Data(), ny, ymin, ymax, nx, xmin, xmax);



  // for (int i=1; i<=nx; ++i) {
  //   for (int j=1; j<=ny; ++j) {
  //     double val = h->GetBinContent(i, j);
  //     double err = h->GetBinError(i, j);
  //     ht->SetBinContent(j, i, val);
  //     ht->SetBinError(j, i, err);
  //   }
  // }////////////++++++++++++++++++++++++++++++++

  ht->SetTitle(h->GetTitle());
  ht->GetXaxis()->SetTitle(ay->GetTitle());
  ht->GetYaxis()->SetTitle(ax->GetTitle());
  return ht;  
}
*/

TH2D* Transpose(TH2* h)
{
  static int ihist = 1;
  TString name = Form("%sTranspose%d", h->GetName(), ihist); ihist++;

  TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis(); 
  double xmin = ax->GetXmin(), xmax = ax->GetXmax(); 
  double ymin = ay->GetXmin(), ymax = ay->GetXmax(); 
  int nx = h->GetNbinsX(), ny = h->GetNbinsY();

  // TMatrixD* Aptr = RooUnfoldResponse::H2M(h, nx, ny, 0, kFALSE);
  // TCanvas* ccc = new TCanvas("ccc","ccc",1);
  // Aptr->Draw("colz");


  TH2D* ht = new TH2D(name.Data(), name.Data(), ny, ymin, ymax, nx, xmin, xmax);
  for (int i=1; i<=nx; ++i) {
    for (int j=1; j<=ny; ++j) {
      double val = h->GetBinContent(i, j);
      double err = h->GetBinError(i, j);
      ht->SetBinContent(j, i, val);
      ht->SetBinError(j, i, err);
    }
  }
  ht->SetTitle(h->GetTitle());
  ht->GetXaxis()->SetTitle(ay->GetTitle());
  ht->GetYaxis()->SetTitle(ax->GetTitle());
  return ht;  
}

TH1D* InnerProduct(TH2* hR, TH1* hT, TH1* hResult)
{
  // Calculate h[j] = sum_i(hR[i,j] hT[i]) 
  //
  // R must have the same number of x bins as T.
  // R must have the same number of y bins as h.

  TAxis *ax = hR->GetXaxis(), *ay = hR->GetYaxis(); 
  //  TArrayD* ybins = ay->GetXbins(); // somehow wrong, range is 0..1
  double xmin = ax->GetXmin(), xmax = ax->GetXmax(); 
  double ymin = ay->GetXmin(), ymax = ay->GetXmax(); 
  int nx = hR->GetNbinsX(), ny = hR->GetNbinsY();
  int nxT = hT->GetNbinsX();

  TH1D* h = dynamic_cast<TH1D*>(hResult);
  if (h) { // Use provided histogram
    h->Reset();
  }
  else { // Or make a new one, if nothing provided
    static int ihist = 1;
    TString name = Form("hFold%d", ihist); ihist++;
    //  TH1D* h = new TH1D(name.Data(), name.Data(), nx, xbins->GetArray()); // no workie
    h = new TH1D(name.Data(), name.Data(), ny, ymin, ymax);
  }
  
  // Check binning and ranges
  if (nx!=nxT || nx==0 || nxT==0) { // hT must have nx bins
    gROOT->Error("InnerProduct()", "%d x bins in %s != %d in %s", 
	  nx, hR->GetName(), nxT, hT->GetName());
    return 0;
  }
  if (ny != h->GetNbinsX()) { // Result must have ny bins
    gROOT->Error("InnerProduct()", "%d y bins in %s != %d in %s",
	  ny, hR->GetName(), h->GetNbinsX(), h->GetName());
    return 0;
  }
  if (xmin != h->GetXaxis()->GetXmin()) {
    gROOT->Warning("InnerProduct()", "Matrix xmin %3.2g != result %3.2g",
	    xmin,  h->GetXaxis()->GetXmin());
  }
  if (xmax != h->GetXaxis()->GetXmax()) {
    gROOT->Warning("InnerProduct()", "Matrix xmax %3.2g != result %3.2g",
	    xmax,  h->GetXaxis()->GetXmax());
  }

  // Calculate h[j] = sum_i R[i][j]*T[i]
  for (int j=1; j<=ny; ++j) {
    double val = 0., err = 0.;
    for (int i=1; i<=nx; ++i) {
      val += hR->GetBinContent(i,j) * hT->GetBinContent(i);
      err = 0.;
    }
    h->SetBinContent(j, val);
    //    h->SetBinError(j, err);
    h->SetBinError(j, TMath::Sqrt(val));
  }


  return h;
}

TH2D* MatrixProduct(TH2* hA, TH2* hB, TH2* hResult)
{
  // Calculate AB[i,j] = sum_k(A[i,k] B[k,j]) 
  //
  // Required:
  //
  // nxa = nyb  (# cols in A = rows in B)
  // nya = nyab (# rows in A = rows in AB)
  // nxb = nxab (# cols in B = cols in AB)

  //  TAxis *ax = hA->GetXaxis(); 
  TAxis* ay = hA->GetYaxis(); 
  TAxis *bx = hB->GetXaxis();//, *by = hB->GetYaxis(); 
  double 
    //    axmin = ax->GetXmin(), axmax = ax->GetXmax(),
    bxmin = bx->GetXmin(), bxmax = bx->GetXmax(), 
    aymin = ay->GetXmin(), aymax = ay->GetXmax();
    //bymin = by->GetXmin(), bymax = by->GetXmax(); 
  int 
    nxa = hA->GetNbinsX(), nya = hA->GetNbinsY(),
    nxb = hB->GetNbinsX(), nyb = hB->GetNbinsY();

  TH2D* h = dynamic_cast<TH2D*>(hResult);
  if (h) { // Use provided histogram
    h->Reset();
  }
  else { // Or make a new one, if nothing provided
    static int ihist = 1;
    TString name = Form("hProd%d", ihist); ihist++;
    h = new TH2D(name.Data(), name.Data(), nxb, bxmin, bxmax, nya, aymin, aymax);
  }
  
  // Check binning and ranges
  if (nxa!=nyb || nxa==0 || nyb==0) {
    gROOT->Error("MatrixProduct()", "%d x bins in %s != %d y bins in %s", 
	  nxa, hA->GetName(), nyb, hB->GetName());
    return 0;
  }
  if (nxb != h->GetNbinsX()) {
    gROOT->Error("MatrixProduct()", "%d x bins in %s != %d in %s",
	  nxb, hB->GetName(), h->GetNbinsX(), h->GetName());
    return 0;
  }
  if (nya != h->GetNbinsY()) {
    gROOT->Error("MatrixProduct()", "%d y bins in %s != %d in %s",
	  nya, hA->GetName(), h->GetNbinsY(), h->GetName());
    return 0;
  }
  if (bxmin != h->GetXaxis()->GetXmin()) {
    gROOT->Warning("MatrixProduct()", "Matrix xmin %3.2g != result %3.2g",
	    bxmin,  h->GetXaxis()->GetXmin());
  }
  if (bxmax != h->GetXaxis()->GetXmax()) {
    gROOT->Warning("MatrixProduct()", "Matrix xmax %3.2g != result %3.2g",
	    bxmax,  h->GetXaxis()->GetXmax());
  }

  // Calculate AB[i,j] = sum_k(A[i,k] B[k,j]) 
  for (int i=1; i<=nxb; ++i) {
    for (int j=1; j<=nya; ++j) {
      double val = 0., err = 0.;
      for (int k=1; k<=nxa; ++k) {
	val += hA->GetBinContent(i,k) * hB->GetBinContent(k,j);
	err = 0.;
      }
      h->SetBinContent(i,j,val);
      h->SetBinError(i,j,err);
    }
  }
  return h;
}

double UnfoldingScaleFactor(TH2* h, int definition)
{
  // Return normalization factor used to guarantee convergence. This
  // is different than K_{\eta,\rho} in arxiv:1111.3387, I found it in
  // numerical recipes 3rd ed. section 2.5.1 (iterative improvement)
  int nx = h->GetNbinsX(), ny = h->GetNbinsY();

  if (definition==0) {
    double max_i=0, max_j=0;
    // Set max_i
    for (int i=1; i<=nx; ++i) {
      double sum = 0;
      for (int j=1; j<=ny; ++j) {
	sum += TMath::Abs(h->GetBinContent(i, j));
      }
      if (max_i < sum) max_i = sum;
    }
    // Set max_j
    for (int j=1; j<=ny; ++j) {
      double sum = 0;
      for (int i=1; i<=nx; ++i) {
	sum += TMath::Abs(h->GetBinContent(i, j));
      }
      if (max_j < sum) max_j = sum;
    }
    return max_i * max_j;
  }
  
  // This is usually a smaller factor. It leads to smaller differences
  // between steps, and slower convergence. Sometimes this is desired
  // if a very small number of unfolding iterations is being used.
  if (definition==1) {
    double sum = 0;
    for (int i=1; i<=nx; ++i) {
      for (int j=1; j<=ny; ++j) {
	double aij = h->GetBinContent(i, j);
	sum += aij*aij;
      }
    }
    return sum;
  }

  if (definition==2) {
    // double sum = 0;
    // for (int i=1; i<=nx; ++i) {
    //   for (int j=1; j<=ny; ++j) {
    // 	double aij = h->GetBinContent(i, j);
    // 	sum += aij*aij;
    //   }
    // }
    // return sum;
  }
  

  return 0.0;
}


void ArrayFromTH2(const TH2* h, double val[], double err[])
{
  int nx = h->GetNbinsX(), ny = h->GetNbinsY();
  int ip = 0;
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      val[ip] = h->GetBinContent(i+1, ny-j);
      err[ip] = h->GetBinError  (i+1, ny-j); 
      ip++;
    }
  }
  return;
}

void ArrayFromTH1(const TH1* h, double val[], double err[])
{
  int nx = h->GetNbinsX();
  int ip = 0;
  for (int i=0; i<nx; i++) {
    val[ip] = h->GetBinContent(i+1);
    err[ip] = h->GetBinError  (i+1); 
    ip++;
  }
  return;
}

void MatrixFromTH2(const TH2& h, TMatrixD& A)
{
  const int nx = h.GetNbinsX(), ny = h.GetNbinsY();
  const int nbins = nx*ny;
  double a[nbins];
  int ip = 0;

  A.ResizeTo(ny, nx);

  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      a[ip] = h.GetBinContent(i+1, ny-j);
      ip++;
    }
  }
  A.SetMatrixArray(a);
}

void VectorFromTH1(const TH1& h, TVectorD& v)
{
  int nb = h.GetNbinsX();
  v.ResizeTo(nb);
  for (Int_t i= 0; i < nb; i++) {
    v(i) = h.GetBinContent(i+1);
  }
}

void HistFromVector(const TVectorD& v, TH1* h)
{
  h->Reset();
  int nb = h->GetNbinsX();

  for (int i = 0; i < nb; i++) {
    h->SetBinContent(i+1, v(i));
  }
}

//==============================================================================
// Chi squared
//==============================================================================

Double_t ChiSquared(TH1* h, TFormula hyp, TString opt = "")
{
  // Calculate chi2 wrt hypothesis fn.
  double chi2 = 0;
  for (int j=1; j<=h->GetNbinsX(); j++) {
    double xj = h->GetBinContent(j)*h->GetBinWidth(j);
    double xt = hyp.Eval(h->GetBinCenter(j));
    double err = h->GetBinError(j);
    if (err==0 && xj!=0) {
      gROOT->Warning("ChiSquared()", 
		     "Bin %d content = %3.2g but error = 0", j, xj);
    }

    Printf("xj, xt, err = %3.2g, %3.2g, %3.2g", xj, xt, err);
    chi2 += (err!=0) ? (xj-xt)*(xj-xt)/err/err : 0.;
  }

  // If you want total, not chi2/NDF
  if (opt.Contains("TOT") || opt.Contains("tot") )
    return chi2;
  
  return chi2 / h->GetNbinsX(); // TODO should be N-1? look up....
}

Double_t ChiSquared(TH1* h, TH1* hTrue, TString opt = "")
{
  // Calculate chi2 wrt hTrue (assumed to have no error).
  double chi2 = 0;
  int nbh = h->GetNbinsX();
  int nbt = hTrue->GetNbinsX();
  if (nbh != nbt) {
    gROOT->Error("ChiSquared()", 
		 "bin mismatch %d in h, %d in hTrue", nbh, nbt);
    return -1.;
  }

  for (int j=1; j<=nbh; j++) {
    double xj = h->GetBinContent(j);
    double xt = hTrue->GetBinContent(j);
    double err = h->GetBinError(j);
    if (err==0 && xj!=0) {
      gROOT->Warning("ChiSquared()", 
		     "Bin %d content = %3.2g but error = 0", j, xj);
    }

    if (0)
      Printf("xj, xt, err = %3.2g, %3.2g, %3.2g", xj, xt, err);
    chi2 += (err!=0) ? (xj-xt)*(xj-xt)/err/err : 0.;
  }

  // If you want total, not chi2/NDF
  if (opt.Contains("TOT") || opt.Contains("tot") )
    return chi2;
  
  return chi2 / h->GetNbinsX(); // TODO should be N-1? look up....
}

TGraph* Chi2vsStep(TObjArray* hists, TH1* hTrue)
{
  //  Calculate chi2 wrt hTrue vs sampled step
  TGraph* g = new TGraph();

  for (int i=0; i<hists->GetEntries(); i++) {
    TH1* h = (TH1*)hists->At(i);
    if (h) {
      // double chi2 = ChiSquared(h, hTrue);
      //      cout << chi2 << endl;
      g->SetPoint(i, i, ChiSquared(h, hTrue));
    }
  }
  
  // For looks
  g->SetLineColor(kRed);
  g->SetMarkerColor(kBlack);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(0.75);
  
  //    PlotUtils::set_tgraph_props(g, kRed, kBlack, kFullCircle, 0.75);
  
  return g;
}

void ScaleToValue(TH1* h, double norm, TString opt)
{
  h->Scale(norm/h->Integral());
  if (opt.Contains("e")) { // sqrt(n) uncertainty
    for (int j=1; j<=h->GetNbinsX(); j++) {
      double val = h->GetBinContent(j);
      double err = (val>=1) ? sqrt(val) : 1.;
      h->SetBinError(j, err);
    }  
  }
}

TH1* MinL2NormHist(TObjArray* xhists, TGraph* g)
{
  // Return histo from xhists with minimum norm, according to g created by L2NormVsStep().

  TH1* h = 0;
  double y = 0, ymin = 1234567890;
  double kmin = 0;
  for (int i=0; i<g->GetN(); i++) {
    y = g->GetY()[i];
    if (y < ymin) {
      ymin = y;
      kmin = i;
    }
  }
  h = (TH1*)xhists->At(kmin);
  if (h) 
    return h;
  
  gROOT->Warning("MinL2NormHist()", "Returning empty histogram");
  return h;
}

TH1* BestLCurveHist(TObjArray* xhists, TGraph* g)
{
  // Return histo from xhists with minimum distance from zero on
  // L-Curve.

  TH1* h = 0;
  double x, y;
  double r2 = 0, rmin = 1e99;
  int kmin = 0;
  for (int i=1; i<g->GetN(); i++) {
    x = g->GetX()[i];
    y = g->GetY()[i];
    r2 = x*x + y*y;
    if (r2 < rmin) {
      rmin = r2;
      kmin = i;
    }
  }
  Printf("Best k = %d", kmin);
  h = (TH1*)xhists->At(kmin);
  if (h) 
    return h;
  
  gROOT->Warning("MinL2NormHist()", "Returning empty histogram");
  return h;
}

TGraph* L2NormVsStep(TH1* hTrue, TObjArray* xhists, TString /*opt*/)
{
  // 2-norm ratio ||xn-xtrue|| / ||xtrue|| vs iteration step.
  // Expecting hResp in true-y-vs-meas-x RooUnfold format.

  if (!hTrue) {
    gROOT->Error("L2NormVsStep()", "hTrue is a null pointer");
    return 0;
  }

  TGraph* g = new TGraph();
  int nx =  hTrue->GetNbinsX();
  TVectorD xt = Hist2Vec(hTrue, nx, kFALSE);
  double xtnorm = TMath::Sqrt(xt*xt);
  if (xtnorm<=0) {
    gROOT->Error("L2NormVsStep()", "xtnorm = %3.2g", xtnorm);
    return 0;
  }

  double yval = 0;  
  for (int i=0; i<xhists->GetEntries(); i++) {
    TH1* h = (TH1*)xhists->At(i);

    if (h) {
      if (h->GetNbinsX() != nx) {
	gROOT->Error("L2NormVsStep()",
		     "%d bins in true != %d bins in %s",
		     nx, h->GetNbinsX(), h->GetName());
	return 0;
      }
      TVectorD xn = Hist2Vec(h, nx, kFALSE);
      TVectorD r(xn);
      r -= xt;
      yval = TMath::Sqrt(r.Norm2Sqr())/xtnorm;
      g->SetPoint(i, i, yval);
    }
  }
  
  // For looks
  g->SetLineColor(kBlack);
  g->SetMarkerColor(kBlack);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1.5);
  g->SetLineWidth(2);
  g->SetTitle("Error history;iteration step k;"
	      "||x^{(k)}-x^{(true)}||_{2}/||x^{(true)}||_{2}"); 
  return g;
}

TMatrixD Toeplitz(int m1, int n1, double col[], double row[])
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
    gROOT->Warning("Toeplitz()", 
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


TMatrixD* GetDerivativeMatrix(int n, int d, TString opt)
{
  // Return a d-th derivative operator matrix.
  // Default size is n-d x n.
  // Unit matrix is returned if d=0.

  int nd = n-d;

  TMatrixD* L = new TMatrixD(nd, n);
  TVectorD c(d+1);
  c.Zero();
  
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

  //  c.Print();

  if (opt.Contains("bczero") || opt.Contains("bcrefl")) {
    L->ResizeTo(n,n);
    int k;
    for (int i=0; i<n; i++) {
      for (int j=0; j<d+1; j++) {
	k = j+i-1;
	if (k<0 || k>=n)
	  continue;
	(*L)(i,k) = c(j);
      }
    }

    (*L)(n-1, n-1) = -1;
    
    if (opt.Contains("bcrefl")) {
      // reflexive bc - 
      // just change UL and LR corners -2 --> -1
      (*L)(0,0) = -1;
      (*L)(n-1, n-1) = -1;

      TDecompQRH decomp(*L);
      TMatrixD R = decomp.GetR();
      for (int i=0; i<nd; i++) {
	for (int j=0; j<d+1; j++) {
	  (*L)(i,j) = R(i,j);
	}
      }
    }
  }
  else {
    for (int i=0; i<nd; i++) {
      for (int j=0; j<d+1; j++) {
	(*L)(i,j+i) = c(j);
      }
    }
  }
  
  return L;
}

TMatrixD* GetNullSpaceMatrix(int n, int d)
{
  TMatrixD* W = new TMatrixD(n, d);

  // Compute basis of n x d null space via
  // modified Graham-Schmidt orthogonalization
  for (int i=0; i<n; i++)
    for (int j=0; j<d; j++)
      (*W)(i, j) = j ? TMath::Power(i+1, j) : 1;
  
  for (int j=0; j<d; j++) {
    TVectorD colj(TMatrixDColumn((*W),j));

    for (int i=0; i<j; i++) {
      TVectorD coli(TMatrixDColumn((*W),i));
      colj -= (colj*coli) / (coli*coli) * coli;
    }
    colj *= 1./TMath::Sqrt(colj.Norm2Sqr());

    for (int i=0; i<n; i++) {
      (*W)(i,j) = colj(i);    
    }
  }

  return W;

}

TMatrixD Null(TMatrixD& A)
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

void LTSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y)
{
  // Computes  x = (L_p)'*y
  // where L_p is the generalized inverse of L (aka L_A^dagger)
  
  int p = L.GetNrows();
  int n = L.GetNcols();
  int ly = y.GetNoElements();
  
  // if (result.GetNoElements() != p)  
  //   result.ResizeTo(p);
  
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

    // update y to new size ?
    // y = y1;

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


void LSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y, 
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

    //    result.ResizeTo(x.GetNoElements());
    result = x;
    return;
  }
  if (p > n) {
    TMatrixD Ltmp = L;
    TMatrixD Linv = MoorePenroseInverse(Ltmp); // now Linv is n x p (?)
    TVectorD x = Linv * y;
    //    result.ResizeTo(x.GetNoElements());
    result = x;
    return;
  }
  else {
    // x = L(:,1:p)\y;
    // x = [x;zeros(nu,ly)] - W*(T(:,1:p)*x);
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

// void LSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y, 
// 	    const TMatrixD& W, const TMatrixD& T)
// {
//   // Computes  x = L_p*y
//   // where L_p is the A-weighted generalized inverse of L.
//   int p = L.GetNrows();
//   int n = L.GetNcols();
//   int nu = n-p;
//   //  int ly = y.GetNoElements();

//   if (result.GetNoElements() != n)  
//     result.ResizeTo(n);

//   result.Zero();

//   if (nu == 0) { // Square L
//     // x = L\y
//     TMatrixD Linv(L);
//     // Linv.Invert();
//     // result = Linv*y;
//     TDecompSVD decomp(L);
//     bool ok = false;
//     result = decomp.Solve(y, ok);
//     if (!ok)
//       gROOT->Warning("LSolve()", "Failed to solve L*x = y");
//     return;
//   }
//   else {
//     // x = L(:,1:p)\y;
//     // x = [x;zeros(nu,ly)] - W*(T(:,1:p)*x);

//     TVectorD ysub(y.GetSub(0, p-1));
//     TMatrixD Lsub(L.GetSub(0, p-1, 0, p-1));
//     TMatrixD Tsub(T.GetSub(0, nu-1, 0, p-1));

//     // Printf("ysub (%d)", ysub.GetNrows());
//     // Printf("Lsub (%d x %d)", Lsub.GetNrows(), Lsub.GetNcols());
//     // Printf("T (%d x %d)", T.GetNrows(), T.GetNcols());

//     TDecompSVD decomp(Lsub);
//     bool ok = false;
//     TVectorD xsub = decomp.Solve(ysub, ok);
//     if (!ok)
//       gROOT->Warning("LSolve()", "Failed to solve L*x = y");

//     // Lsub.Invert();
//     // TVectorD xsub = Lsub*ysub;
//     //    Printf("xsub (%d), nu=%d", xsub.GetNoElements(), nu);

//     result.SetSub(0, xsub);
//     //    result.Print();
//     //    result -= W*Tsub*xsub;
//     result -= W*(Tsub*xsub);
//   }
//   return;
//  }

TMatrixD MoorePenroseInverse(TMatrixD& A, double tol)
{
  // Compute pseudoinverse of A as V*diag(1/sigma_i)*U' (Numerical Recipes eq 2.6.6)
 int m = A.GetNrows(), n = A.GetNcols();
 int max = m;
 if (n>max) 
   max = n;
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

 TMatrixD Ainv = V*Siginv*UT;
 Ainv.ResizeTo(n,m);

 return Ainv;
}


TMatrixD* GetPseudoInverse(TMatrixD& A)
{
  int m = A.GetNrows(), n = A.GetNcols();
  TMatrixD Amm(m,m);
  Amm.SetSub(0, 0, A);
  
  TDecompSVD decomp(Amm);
  double c = decomp.Condition();
  gROOT->Info("GetPseudoInverse()", "condition = %3.2g", c); 
  bool status = decomp.Invert(Amm); 
  gROOT->Info("GetPseudoInverse()", "status = %d", (int)status); 

  // Now the square matrix Smm has been pseudo-inverted, but contains
  // many empty rows.

  TMatrixD* Ainv = new TMatrixD(n,m);
   for (int i=0; i<n; i++) {
     TMatrixDRow(*Ainv, i) = TMatrixDRow(Amm, i);
   }
   return Ainv;
}

// TGraphTime* Animation(TObjArray* moveObjs, TObjArray* statObjs, TString opt)
// {
//   static int iAnim = 0; iAnim++;
//   int nFrames = moveObjs->GetEntries();
//   Info("Animation()", "Creating %d frame sequence...", nFrames);
//   TGraphTime* anim = new TGraphTime(nFrames,0,0,1,1);
//   anim->SetName(Form("anim%d", iAnim));

//   for (int n=0; n<nFrames; n++) {

//     // Add stationary objects to this frame
//     if (statObjs) {
//       for (int i=0; i<statObjs->GetEntries(); i++) {
// 	TObject* sta = statObjs->At(i);
// 	TString drawOpt = i ? "same" : "";
// 	if (!opt.IsNull())
// 	  drawOpt += opt;
// 	anim->Add(sta, n, drawOpt);
//       }
//     }
//     // Add changing objects
//     TObject* mov = moveObjs->At(n);
//     TString drawOpt2 = statObjs ? "same" : "";
//     if (!opt.IsNull())
//       drawOpt2 += opt;
//     anim->Add(mov, n, drawOpt2);
//   }
//   anim->SetTitle("animation");
//   anim->SetSleepTime(100); // ms (default = 0)
//   return anim;
// }

TH1* RatioHist(TH1* h1, TH1* h2)
{
  TH1* h = (TH1*)h1->Clone(Form("%s_ratio", h1->GetName()));
  h->Divide(h2);
  return h;
}

TObjArray* RatioHists(TObjArray* hists, TH1* denom)
{
  TObjArray* ratios = new TObjArray();
  for (int i=0; i<hists->GetEntries(); i++) {
    TH1* h = (TH1*)hists->At(i);

    if (h && h->InheritsFrom("TH1")) {
      TH1* hr = (TH1*)h->Clone(Form("%s_ratio", h->GetName()));
      hr->Divide(denom);
      //      cout << hr->GetName() << endl;
      ratios->Add(hr);
    }
    else
      gROOT->Warning("RatioHists()", "Null or non-TH1 pointer");
  }
  
  return ratios;
}

void ToyJetSpectraModel(const int m, const int n, 
			const double x1, const double x2, 
			TMatrixD& A, TVectorD& x, TVectorD& b, 
			double nCounts)
{
  // A must have size m x n
  // x must have length n
  // b must have length m

  if (n%2) {
    gROOT->Error("Shaw()", "Even binning required");
    return;
  }

  TMatrixD Adet(m,n);
  TMatrixD Abkg(m,n);

  // For computing A
  double h = (x2-x1)/n;
  
  // For computing x vector
  double npwr = 4.0, T = 1, ptj = 0.;
  double plaw[n];
  
  for (int j=0; j<n; j++) {
    ptj  = x1 + (j + 0.5)*h;    
    // x is "generated" or "truth" vector
    plaw[j] = ptj >= 0 ? ptj*TMath::Power(1+ptj/npwr/T, -npwr) : 0.0;
    x(j) = plaw[j];
  }

  // Normalize to nCounts
  x = nCounts/x.Sum() * x;

  double apar = 0.0, bpar = 0.1, cpar = 0.01,
    mean=0, sigma=0, pti=0;

  double sigma_bkg = 10;  

  for (int j=0; j<n; j++) { // col index
    mean = x1 + (j + 0.5)*h;
    for (int i=0; i<m; i++) { // row index

      // true pt at point i
      pti  = x1 + (i + 0.5)*h;

      // sigma of detector response
      sigma = pti > 0.0 ? apar + sqrt(bpar*fabs(mean)) + cpar*fabs(mean) : apar;
      // double res = plaw[j]*TMath::Gaus(pti, mean, sigma, true);
      double res = 
	pti*TMath::Power(1+pti/4/6, -4)
	*TMath::Landau(pti, mean, sigma, true);

      // background fluctuations
      double bkg = TMath::Gaus(pti, mean, sigma_bkg, true);
      //      A(i,j) = val > 1e-15 ? val : 0.;
      Adet(i,j) = res;
      Abkg(i,j) = bkg;
    }
  }

  A = Abkg*Adet;

  // Normalize A and calculate b
  b = A*x;
  A *= nCounts/b.Sum();
  b = A*x;

  return;
}

void ShawSystem(const int n, TMatrixD& A, TVectorD& x, TVectorD& b)
{

  // Create a complete test problem for unfolding exercises that is
  // mathematically identical to shaw.m in the Matlab "regularization
  // tools" examples.
  //
  // Just do, e.g.
  // TMatrixD A(n,n);
  // TVectorD xt(n), bt(n);      // true x and b (no noise)
  // ShawSystem(n, A, xt, bt);
  // Where is an even integer.

  if (n%2) {
    gROOT->Error("Shaw()", "Even binning required");
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

  return;
}
/*
TH1* UnfoldTSVD(const TMatrixD& A, 
		const TVectorD& b, 
		const int nsv, 
		TObjArray* hists, 
		TObjArray* extras, 
		const double x1, 
		const double x2,
		int ireturn=nsv-1)
{
  TDecompSVD decomp(A);
  TVectorD sigma_vec = decomp.GetSig();
  TMatrixD UT = decomp.GetU(); UT.T();
  TMatrixD V = decomp.GetV();
  TVectorD x(sigma_vec.GetNoElements());

  for (int i=0; i<nsv; i++) {
    TVectorD uTi = TMatrixDColumn(UT, i);
    double sig = sigma_vec(i);
    double svd_coeff = sig ? (uTi*b)/sig : 0.;
    TVectorD v =  TMatrixDColumn(V, i);
    x += svd_coeff*v;

    if (hists) // Save ith partial sum as a TH1D
      hists->Add(XHist(x,Form("PCGLS%d",ihist),i,x1,x2,normTo,""));
  }
  return hists->At(ireturn);
}
*/
void SVDAnalysis(TMatrixD& A, TVectorD& b, TObjArray* output)
{

  // Decompose A as U*Sigma*V' and study the components. If A is m x
  // n, then U is a column-orthogonal m x n matrix, Sigma is a
  // diagonal n x n matrix (represented below as a vector), and V is
  // a column-orthonormal n x n matrix (V'*V = 1).

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

  TH1D *hSigma=0, *hUb=0, *hUbSig=0;

  hSigma     = new TH1D("hSigma", "#sigma_{i} ", nvals, 0, nvals);
  hUb    = new TH1D("hUb", "|u^{T}_{i}*b| ", mU, 0, mU);
  hUbSig = new TH1D("hUbSig", "|u^{T}_{i}*b| / #sigma_{i} ", mU, 0, mU);
  SetHistProps(hSigma, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(hUb, kBlue, kNone, kBlue, kFullSquare, 1.0);
  SetHistProps(hUbSig, kRed, kNone, kRed, kOpenSquare, 1.0);
  //  RooUnfoldResponse::V2H(sigma_vec, hSigma, nvals);
  Vec2Hist(sigma_vec, hSigma, nvals);
  
  output->Add(hSigma);
  output->Add(hUb);
  output->Add(hUbSig);

  for (int i=0; i<nU; i++) {
    TVectorD ui = TMatrixDColumn(U, i);
    TH1D* hui = new TH1D(Form("u_%d", i), Form("u_{%d}", i), mU, 0, mU);
    //    RooUnfoldResponse::V2H(ui, hui, mU);
    Vec2Hist(ui, hui, mU);
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

void SetHistProps(TH1* h,
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

TMatrixD LMatrix(const int n, const int kind)
{

  // Return p x n smoothing matrix. The number of rows p is
  // case-dependent.

  // Top row and left column.
  // Make sure that row[0] = col[0]
  double row[n];
  for (int j=0; j<n; j++) row[j]=0.0;

  if (kind == 0) { // n x n identity matrix. p = n. W=0.
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    return Toeplitz(n,n,col,row);
  }
  if (kind == 1) { // 1st deriv, no BC assumptions. p = n-1. W=const.
    double col[n-1];
    for (int i=0; i<n-1; i++) col[i]=0.0;
    row[0] = col[0] = -1;
    row[1] = 1;
    return Toeplitz(n-1,n,col,row);
  }
  if (kind == 2) { // 2nd deriv, no BC assumptions. p = n-2. W=const, linear.
    double col[n-2];
    for (int i=0; i<n-2; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    row[1] = -2;
    row[2] = 1;
    return Toeplitz(n-2,n,col,row);
  }
  if (kind == 3) { // 1st deriv, BC=0 on L and R. p = n+1. W=0.
    double col[n+1];
    for (int i=0; i<n+1; i++) col[i]=0.0;
    row[0] = col[0] = 1;
    col[1] = -1;
    return Toeplitz(n+1,n,col,row);
  }
  if (kind == 4) { // 2nd deriv, BC=0 on L and R. p = n. W=0.
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = -2;
    row[1] = col[1] = 1;
    return Toeplitz(n,n,col,row);
  }
  if (kind == 5) { // 1st deriv, reflect at L,R. p = n-1. W=const. Same as 1.
    double col[n-1];
    for (int i=0; i<n-1; i++) col[i]=0.0;
    row[0] = col[0] = -1;
    row[1] = 1;
    return Toeplitz(n-1,n,col,row);
  }
  if (kind == 6) { // 2nd deriv, reflect at L,R. p = n. W=const.
    double col[n];
    for (int i=0; i<n; i++) col[i]=0.0;
    row[0] = col[0] = -2;
    row[1] = col[1] = 1;
    TMatrixD L = Toeplitz(n,n,col,row);
    L(0,0) = -1;
    L(n-1,n-1) = -1;
    return L;
  }

  TMatrixD Ldefault(n,n);
  Ldefault.UnitMatrix();
  return Ldefault;
}

TH2D* M2H(TMatrixD& A, TString hName, double x1, double x2, double y1, double y2)
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

TMatrixD Hist2Matrix(const TH2* h)
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


TMatrixD Hist2Matrix(const TH2* h, Int_t nx, Int_t ny, const TH1* norm, Bool_t overflow)
{
  // Returns Matrix of values of bins in a 2D input histogram
  Int_t first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD m(nx, ny);
  if (!h) return m;
  for (Int_t j= 0; j < ny; j++) {
    Double_t nTrue= norm ? GetBinContent(norm, j, overflow) : 1.0;
    if (nTrue == 0.0) {
      for (Int_t i= 0; i < nx; i++) {
        m(i,j)= 0.0;
      }
    } else {
      for (Int_t i= 0; i < nx; i++) {
        m(i,j)= h->GetBinContent(i+first,j+first) / nTrue;
      }
    }
  }
  return m;
}

void Vec2Hist(const TVectorD& v, TH1* h, Int_t nb, Bool_t overflow)
{
  // Sets the bin content of the histogram as that element of the input vector
  h->Reset();  // in particular, ensure under/overflows are reset
  if (overflow) nb += 2;
  for (Int_t i= 0; i < nb; i++) {
    Int_t j= GetBin(h, i, overflow);
    h->SetBinContent(j, v(i));
  }
}

TVectorD Hist2Vec(const TH1* h, Int_t nb, Bool_t overflow)
{
  // Returns TVectorD of the bin contents of the input histogram
  if (overflow) nb += 2;
  TVectorD v(nb);
  if (!h) return v;
  for (Int_t i= 0; i < nb; i++) {
    v(i) = GetBinContent(h, i, overflow);
  }
  return v;
}

double GetBinContent(const TH1* h, Int_t i, Bool_t overflow)
{
  // Bin content by vector index
  return h->GetBinContent(GetBin(h, i, overflow));
}

int GetBin(const TH1* h, Int_t i, Bool_t overflow)
{
  // vector index (0..nx*ny-1) -> multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
  return (h->GetDimension()<2) ? i+(overflow ? 0 : 1) : GetBinDim(h,i);
}

int GetBinDim(const TH1* h, Int_t i)
{
  // Copied from RooUnfoldResponse.  Converts from vector index
  // (0..nx*ny-1) or (0..nx*ny*nz-1) to multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) or
  // (0..(nx+2)*(ny+2)*(nz+2)-1), skipping under/overflow bins.
  Int_t ndim= h->GetDimension(), nx= h->GetNbinsX();
  if        (ndim == 2) {
    return (i%nx+1) + (nx+2)*(i/nx+1);
  } else if (ndim == 3) {
    Int_t ny= h->GetNbinsY();
    return (i%nx+1) + (nx+2)*((i/nx)%ny+1 + (ny+2)*(i/(nx*ny)+1));
  }
  return i+1;   // not used
}

