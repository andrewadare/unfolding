#ifndef UnfoldingUtils_h
#define UnfoldingUtils_h

class TH1;
class TH1D;
class TH2;
class TH2D;
class TH3;
class TH3D;
class TString;
class TGraph;
class TCanvas;
class TF1;

#include "TMatrixD.h"
#include "TVectorD.h"

struct SVDResult         // Output from SVDAnalysis().
{
  TH1D* sigma;           // Singular value spectrum.
  TH1D* UTb;             // Left sing. vectors * b.
  TH1D* UTbAbs;          // |U'*b| (mainly for plotting)
  TH1D* coeff;           // SV coefficients uTb/sigma
  TH2D* U;               // Left sing. vectors
};

class GSVDResult        // Output from GSVDAnalysis().
{
 public:

  GSVDResult();
  virtual ~GSVDResult() {}
  
  int n;                 // Column count of A (or L)
  int m;                 // Row count of A
  int p;                 // Row count of L
  double lambda;         // Regularization parameter used
  TVectorD alpha;        // alpha(i<n-p) = 1, else C(n-p,n-p)..C(n,n).
  TVectorD beta;         // beta(i<n-p)  = 0, else S(n-p,n-p)..S(n,n).
  TVectorD gamma;        // alpha/beta. Valid for gamma(i >= n-p).
  TVectorD f;            // Tik. filter factors. Size n. f(i<n-p) = 1.
  TVectorD UTb;          // Left sing. vectors * b. (n)
  TVectorD coeff;        // GSV coeffs uT*b/alpha. (n)
  TVectorD regc;         // Regularized (filtered) coeffs. (n)
  TMatrixD X;            // Columns = GSVD basis vectors (n x n)
  TMatrixD U;            // Columns = left sing. vectors of A (m x m)
  TMatrixD V;            // Columns = left sing. vectors of L (p x p)
  TMatrixD L;            // Smoothing matrix used (p x n)
  TMatrixD A;            // Coefficient matrix used (m x n)
  TMatrixD Ap;           // A^#: regularized inverse of A (n x m)
  TMatrixD covw;         // Covariance matrix of wreg Ap * Ap' (n x n)
  TMatrixD covx;         // xini * covw * xini (n x n)
  TMatrixD covxInv;      // Hocker eq. 53 (n x n)
  TMatrixD covb;         // Either fMatb or I (m x m)
  TVectorD b;            // Measured RHS vector used (m)
  TVectorD bInc;         // Incompatible b component (I-UU')b (m x m)
  TVectorD breg;         // Refolded solution

  TH2D* UHist;           // Left sing. vectors
  TH2D* XHist;           // GSVD basis vectors

  // Solution for this lambda
  TH1D* wregHist;        // scaled result w^lambda
  TH1D* xregHist;        // xini_j * w^lambda_j (solution)
  TH1D* bregHist;        // Refolded distribution A*xreg

  // Abs. values (for visual analysis)
  TH1D* UTbAbs;          // Vector of |U'_i*b| values
  TH1D* coeffAbs;        // Vector of |U'_i*b|/alpha_i values
  TH1D* regcAbs;         // Regularized (filtered) coeffs.
  
  ClassDef(GSVDResult, 1);
};

struct UnfoldingResult
{
  // Solution matrix (for iterative methods, or reg. parameter scan).
  // Sequence of results stored in a TH2. Results are slices along x.
  // Bin k along y is the kth solution.
  TMatrixD WReg;
  TMatrixD XReg;
  TMatrixD WRegErr;
  TMatrixD XRegErr;

  // Covariance matrices of w and x (for single best solution)
  TMatrixD wCov;
  TMatrixD xCov;

  TH3D* hwCov;
  TH3D* hxCov;

  TH2D* WRegHist;
  TH2D* XRegHist;


  TGraph* LCurve;  // Parametric curve of ||Lx||_2 vs. ||Ax-b||_2.
  TMatrixD F; // filter factors vs. iteration number

  // Summed Tikhonov Filter factors vs. lambda
  // Optimal point is smallest lambda where 
  // sum(f_i) >= # significant GSVD U'b coeffs. 
  TGraph* FilterSum;
  double lambdaStf;
  int kStf; // Optimal point index on curve
  TH1D* hStf;

  // Curvature of L-Curve
  TGraph* LCurvature;
  double lambdaLcv;
  int kLcv;     // Optimal point index on curve
  TH1D* hLcv; // Best unfolding result (according to L-Curve)

  // Generalized cross-validation curve,
  // and best lambda & iteration / index
  TGraph* GcvCurve;
  double lambdaGcv;
  int kGcv;
  TH1D* hGcv; // Best unfolding result (according to GCV)

  // Mean of global correlation coefficients
  TGraph* RhoCurve;
  double lambdaRho;
  int kRho;
  TH1D* hRho; // Best unfolding result (according to GCV)
  
  ClassDef(UnfoldingResult, 1);
};

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
  TH1D* XHist(TVectorD& x, TString base, int k, double xMin, 
	      double xMax, double normto, TString opt="");
  TGraph* LogCurvature(TGraph* g, const TVectorD& tVec, int& kMax); // log curvature of parametric curve g(t). kMax is point index of greatest curvature.
  TMatrixD RegularizedInverseResponse(GSVDResult* gsvd, double lambda); // A^#
  TMatrixD LMatrix(const int n, const int kind, double eps = 1.e-5); // Matrix to define smoothing seminorm
  void SetTH1Props(TH1* h,
		   Int_t linecolor,
		   Int_t fillcolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize); 
  
  // Unfolding methods:
  // Preconditioned Conjugate Gradients for Least Squares
  UnfoldingResult UnfoldPCGLS(const int nIterations, 
			      int lMatrixType             = k2DerivBCR,
			      TString opt                 = "",
			      const GSVDResult* gsvd      = 0,
			      const TH2* hA               = 0, 
			      const TH1* hb               = 0, 
			      const TH1* hXini            = 0);
  UnfoldingResult UnfoldPCGLS(const int nIterations, 
			      TMatrixD& L,
			      TString opt                 = "",
			      const GSVDResult* gsvd      = 0,
			      const TH2* hA               = 0, 
			      const TH1* hb               = 0, 
			      const TH1* hXini            = 0);
  
  TH1D* UnfoldTLS();


  // Richardson-Lucy algorithm
  UnfoldingResult UnfoldRichardsonLucy(const int nIterations);
  
  // Regularized Hocker/Kartvilishveli SVD algorithm
  TH1D* UnfoldSVD(double lambda, 
		  TObjArray* output             = 0, 
		  TString opt                   = "", /*"BC0", "BCR"*/
		  TH2* hA                       = 0, 
		  TH1* hb                       = 0, 
		  TH1* hXini                    = 0);
  
  // Regularized chi squared minimization algorithm
  UnfoldingResult UnfoldChiSqMin(TVectorD& regWts, TString opt = ""); 
  
  // General-form Tikhonov solver based on generalized SVD
  UnfoldingResult UnfoldTikhonovGSVD(GSVDResult* gsvd,
				     TVectorD& lambda, 
				     TString opt = "");
  
  // Support functions for chi squared minimization method
  double SmoothingNorm(TVectorD& x, int regtype);
  double Curvature(TVectorD& x);
  double RegChi2(const double *pars);
  
  // Analysis methods
  TGraph* ResidualNorm(TObjArray* hists, double stepSize = 1.);
  SVDResult SVDAnalysis(TH2* hA=0, TH1* hb=0, TString opt="");
  GSVDResult* GSVDAnalysis(TMatrixD& L, double lambda=0, TH2* hA=0, TH1* hb=0, TString opt="");
  TCanvas* DrawSVDPlot(SVDResult, double ymin, double ymax, TString opt="");
  TCanvas* DrawGSVDPlot(GSVDResult*, double ymin, double ymax, TString opt="");
  
  TH2D* UnfoldCovMatrix(int nTrials, 
			int algo, 
			double regPar, 
			TString opt,
			TObjArray* bHists = 0);
  
  // Set and get methods
  void SetTrueRange(double x1, double x2);
  void SetMeasRange(double x1, double x2);
  void SetVerbosity(int val)    {fVerbosity = val;}  
  void SetRegType(int val)      {fRegType   = val;}
  void SetRegWeight(double val) {fRegWeight = val;}
  void SetLMatrix(TMatrixD& L)  {fMatL.ResizeTo(L); fMatL = L;}
  void SetPrior(TH1D* h);
  void SetSmoothingWeights(TVectorD& wts) {fSmoothingWeight = wts;}

  Double_t GetTrueX1()      const {return fTrueX1;}
  Double_t GetTrueX2()      const {return fTrueX2;}
  Double_t GetMeasX1()      const {return fMeasX1;}
  Double_t GetMeasX2()      const {return fMeasX2;}
  Double_t GetbErrNorm();
  Double_t GetbErrMean();
  Double_t GetbErrRMS();
  Int_t GetM()              const {return fM;}         // Number of measured bins
  Int_t GetNMeasBins()      const {return fM;}         // Number of measured bins
  Int_t GetN()              const {return fN;}         // Number of true/generated bins
  Int_t GetNTrueBins()      const {return fN;}         // Number of true/generated bins
  TMatrixD GetA(TString opt = "");                     // "^", "~", or ""
  TVectorD Getb(TString opt = "");                     // "~" or ""
  TMatrixD GetbCovariance() const {return fMatB;}      // Error matrix of b
  TMatrixD GetbBinv()       const {return fMatBinv;}   // Inverse error matrix of b
  TMatrixD GetLMatrix()     const {return fMatL;}
  TVectorD GetxTrue()       const {return fVecXtrue;}  // b, scaled by its error
  TH2D* GetAProbHist() 	    const {return fHistAProb;}
  TH2D* GetAHist()          const {return fHistA;}
  TH2D* GetATildeHist()     const {return fHistATilde;}
  TH2D* GetBCovHist() 	    const {return fHistMeasCov;}
  TH1D* GetbTildeHist()     const {return fHistbTilde;}
  TH1D* GetXiniHist()       const {return fHistXini;}
  TH1D* GetXTrueHist();
  TH1D* GetEffHist()        const {return fHistEff;}
  TH1D* GetPrior()          const {return fHistPrior;}

  // Smoothing matrix types (& descriptions of W = Null(L))
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
		     kGSVDAlgo,
		     kRichLucyAlgo,
		     kChi2MinAlgo,
		     kPCGLSAlgo};
  TString Algorithm(const int alg); // Pass in UnfoldingAlgo type

 protected:
  void ComputeRescaledSystem();

  int fVerbosity;
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
  TH1D* fHistPrior;     // Implemented for R-L and chi^2
  
  TMatrixD fMatA;       // Transfer matrix
  TMatrixD fMatAhat;    // Probability matrix (column sums = 1.0)
  TMatrixD fMatATilde;  // Covariance-preconditioned matrix
  TMatrixD fMatL;       // Smoothing matrix
  TMatrixD fMatB;       // Data covariance matrix
  TMatrixD fMatBinv;    // Data covariance matrix inverse
  TVectorD fVecb;       // Measurement
  TVectorD fVecbErr;    // Measurement errors
  TVectorD fVecbTilde;  // Meas. points / errors
  TVectorD fVecXtrue;   // True solution (for test problems)
  TVectorD fVecXini;    // Model truth distribution
  TVectorD fSmoothingWeight;

  ClassDef(UnfoldingUtils, 1);
};
#endif
