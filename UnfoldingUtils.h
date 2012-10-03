#ifndef UnfoldingUtils_h
#define UnfoldingUtils_h

class TH1;
class TH1D;
class TH2;
class TH2D;
class TString;
class TGraph;
class TCanvas;

#include "TMatrixD.h"
#include "TVectorD.h"

class UnfoldingUtils 
{
 public:

  // Class administration
  UnfoldingUtils();
  UnfoldingUtils(TH2D* hA, 
		 TH1D* hMeas, 
		 TH2D* hMeasCov = 0, 
		 TH1D* hXini    = 0, 
		 TH1D* hXtrue   = 0, 
		 TH1D* hEff     = 0);
  virtual ~UnfoldingUtils() {}
  bool BinningOk();
  
  // Conversion methods
  TVectorD Hist2Vec(const TH1* h);
  TMatrixD Hist2Matrix(const TH2* h);
  TH1D* Vec2Hist(const TVectorD& v, Double_t x1, Double_t x2, TString name, TString title="");  
  TH2D* Matrix2Hist(TMatrixD& A, TString hName,
		    double x1, double x2, double y1, double y2);
  TH1D* XHist(TVectorD& x, TString base, int k, double xMin, 
	      double xMax, double normto, TString opt="");
  TH2D* BandedDiagonalMatrix(TH1* hDpt, const int nMeas, const int nTrue);
  
  // Utility methods
  void NormalizeXSum(TH2* hA, TH1* hN=0); // Modify hA in-place
  TH2* TH2Product(TH2* hA, TH2* hB, TString name);
  TH2D* TH2Sub(TH2* h, int bx1, int bx2, int by1, int by2, TString name);
  TMatrixD MoorePenroseInverse(TMatrixD& A, double tol = 1e-15);
  TMatrixD Null(TMatrixD& A);
  TMatrixD Toeplitz(int m1, int n1, double col[], double row[]);
  TMatrixD LMatrix(const int n, const int kind, double eps = 1.e-5);
  TMatrixD DerivativeMatrix(int n, int d);
  TVectorD ElemMult(const TVectorD& x, const TVectorD& y);
  TVectorD ElemDiv (const TVectorD& x, const TVectorD& y);
  TMatrixD MultRowsByVector(const TMatrixD& M, const TVectorD& v);
  TMatrixD DivColsByVector(const TMatrixD& M, const TVectorD& v, bool makeZeroIfNaN=true);
  TMatrixD OuterProduct(TVectorD a, TVectorD b); // a*b'
  void LTSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y);
  void LSolve(TVectorD& result, const TMatrixD& L, const TVectorD& y, 
	      const TMatrixD& W, const TMatrixD& T);
  void SetTH1Props(TH1* h,
		   Int_t linecolor,
		   Int_t fillcolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize); 
  
  // Test problems
  // Shaw: sets up discretized n x n kernel A, as
  // well as true and measured vectors x and b (both size n).
  // C. B. Shaw, Jr., "Improvements of the resolution of an instrument
  // by numerical solution of an integral equation", J. Math. Anal. Appl. 37
  // (1972), 83–112.
  TObjArray* ShawSystem(const int n, double noise=0.);
  void ShawSystem(const int n, TMatrixD& A, TVectorD& x, TVectorD& b,
		  double noise=0);

  // Unfolding methods:
  // Preconditioned Conjugate Gradients for Least Squares
  TH1D* UnfoldPCGLS(const int nIterations, 
		    TObjArray* hists            = 0,
		    TObjArray* extras           = 0, 
		    int lMatrixType             = k2DerivBCR,
		    TString opt                 = "",
		    const TH2* hA               = 0, 
		    const TH1* hb               = 0, 
		    const TH1* hXini            = 0);
  
  // Richardson-Lucy algorithm
  TH1D* UnfoldRichardsonLucy(const int nIterations, 
			     TObjArray* hists   = 0,
			     TObjArray* extras  = 0, 
			     TString opt        = "",
			     const TH1* hXStart = 0,
			     const TH2* hA      = 0,
			     const TH1* hb      = 0,
			     const TH1* hXini   = 0);
  
  // Regularized Hocker/Kartvilishveli SVD algorithm
  TH1D* UnfoldSVD(double lambda, 
		  TObjArray* output             = 0, 
		  TString opt                   = "", /*"BC0", "BCR"*/
		  TH2* hA                       = 0, 
		  TH1* hb                       = 0, 
		  TH1* hXini                    = 0);
  
  // Regularized chi squared minimization algorithm
  TH1D* UnfoldChiSqMin(double regWt, 
		       TObjArray* output        = 0, 
		       TString opt              = "",
		       TH1* hXStart             = 0, 
		       TH2* hA                  = 0, 
		       TH1* hb                  = 0, 
		       TH1* hXini               = 0); 
  
  // Support functions for chi squared minimization method
  double SmoothingNorm(TVectorD& x, int regtype);
  double Curvature(TVectorD& x);
  double RegChi2(const double *pars);
  
  // Analysis methods
  TGraph* ResidualNorm(TObjArray* hists, double stepSize = 1.);
  void  SVDAnalysis(TObjArray* output, TH2* hA=0, TH1* hb=0, TString opt="");
  TCanvas* DrawSVDPlot(TObjArray* svdhists, double ymin, double ymax);
  TCanvas* DrawGSVDPlot(TObjArray* svdhists, double ymin, double ymax, TString opt="");

  TH2D* UnfoldCovMatrix(int nTrials, 
			int algo, 
			double regPar, 
			TString opt);
  
  // Set and get methods  
  void SetTrueRange(double x1, double x2);
  void SetMeasRange(double x1, double x2);
  void SetRegType(int val)      {fRegType   = val;}
  void SetRegWeight(double val) {fRegWeight = val;}

  Double_t GetTrueX1()      const {return fTrueX1;}
  Double_t GetTrueX2()      const {return fTrueX2;}
  Double_t GetMeasX1()      const {return fMeasX1;}
  Double_t GetMeasX2()      const {return fMeasX2;}
  Int_t GetM()              const {return fM;}         // Number of measured bins
  Int_t GetNMeasBins()      const {return fM;}         // Number of measured bins
  Int_t GetN()              const {return fN;}         // Number of true/generated bins
  Int_t GetNTrueBins()      const {return fN;}         // Number of true/generated bins
  TMatrixD GetA()           const {return fMatA;}      // "Counts" matrix
  TMatrixD GetAProb()       const {return fMatAhat;}   // Prob matrix - cols sum to 1.0
  TMatrixD GetATilde()      const {return fMatATilde;} // A, scaled using error of b
  TMatrixD GetbCovariance() const {return fMatB;}      // Error matrix of b
  TMatrixD GetbBinv()       const {return fMatBinv;}   // Inverse error matrix of b
  TVectorD Getb()           const {return fVecb;}      // vector of measured data points
  TVectorD GetbTilde()      const {return fVecbTilde;} // b, scaled by its error
  TVectorD GetxTrue()       const {return fVecXtrue;}  // b, scaled by its error
  TH2D* GetAProbHist() 	    const {return fHistAProb;}
  TH2D* GetAHist()          const {return fHistA;}
  TH2D* GetATildeHist()     const {return fHistATilde;}
  TH2D* GetBCovHist() 	    const {return fHistMeasCov;}
  TH1D* GetbTildeHist()     const {return fHistbTilde;}
  TH1D* GetXiniHist()       const {return fHistXini;}
  TH1D* GetXTrueHist()      const {return fHistXtrue;}
  
  // Smoothing matrix types (W is Null(L))
  enum LType{kUnitMatrix, // n   x n, W=0 
	     k1DerivNoBC, // n-1 x n, W=const 
	     k2DerivNoBC, // n-2 x n, W=const, linear
	     k1DerivBC0,  // n+1 x n, W=0
	     k2DerivBC0,  // n   x n, W=0
	     k1DerivBCR,  // n-1 x n, W=const
	     k2DerivBCR}; // n   x n, W=const

  // Regularization types
  enum RegType{kNoReg,    // add zero to chi squared
	       k2Norm,    // L_2 norm ||x||_2 = sqrt(x*x).
	       kTotCurv}; // From 2nd derivative
 
  enum UnfoldingAlgo{kSVDAlgo,
		     kRichLucyAlgo,
		     kChi2MinAlgo,
		     kPCGLSAlgo};

 protected:
  void ComputeRescaledSystem();

  bool fTilde;
  Int_t fM, fN;
  Double_t fMeasX1, fMeasX2, fTrueX1, fTrueX2;
  Double_t fRegWeight;   // For chi^2 method
  Int_t fRegType;        // For chi^2 method
  TH2D* fHistAProb;      // aka Ahat
  TH2D* fHistA;
  TH2D* fHistATilde;
  TH1D* fHistMeas;
  TH2D* fHistMeasCov;
  TH1D* fHistbTilde;
  TH1D* fHistXini;
  TH1D* fHistXtrue;
  TH1D* fHistEff;
  
  TMatrixD fMatA;
  TMatrixD fMatAhat;
  TMatrixD fMatATilde;
  TMatrixD fMatB;
  TMatrixD fMatBinv;
  TVectorD fVecb;
  TVectorD fVecbTilde;
  TVectorD fVecXtrue;
  TVectorD fVecXini;

  ClassDef(UnfoldingUtils, 1);
};
#endif
