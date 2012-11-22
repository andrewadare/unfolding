bool savePDF = 0;
const int n=32;                     // number of points

TH1D *hTrue=0;
TH1D *hMeas=0, *hMeasI=0;       // Measured vector b, w/ & w/o noise
TH2D *hResp=0, *hRespT=0;       // Response matrix & transpose
TH1D* hSVD=0;

TObjArray* cList = new TObjArray();
TObjArray* histsRL = new TObjArray();
TObjArray* extrasRL = new TObjArray();
TObjArray* histsPCGLS = new TObjArray();
TObjArray* extrasPCGLS = new TObjArray();
TObjArray* svdResid = new TObjArray();
TObjArray* genSvdAna = new TObjArray();
TObjArray* svdAnim = new TObjArray();
TObjArray* statObjs = new TObjArray();

TLatex lt;

void ShawExample()
{
  gStyle->SetOptStat(0);
  lt.SetNDC();
  LoadLibs();

  // Setup (this instance is not the one used for unfolding).
  UnfoldingUtils shaw;
  shaw.SetTrueRange(0.,1); 
  shaw.SetMeasRange(0.,1);

  // Set the level of Gaussian white noise on b (absolute).  Since
  // cov(err on b) = deltab * I, the prewhitening option "~" is
  // unnecessary for this problem (there is no x^ini).  Only SVD
  // components above this level should contribute to a regularized
  // solution.
  double deltab = 0.05;

  TestProblem p = shaw.ShawSystem(n,deltab);
  hResp  = p.Response;
  hTrue  = p.xTruth;
  hMeas  = p.bNoisy;
  hMeasI = p.bIdeal;

  statObjs->Add(hMeasI); 
  statObjs->Add(hMeas); 
  statObjs->Add(hTrue);
  SetHistProps(hTrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(hMeas, kBlue, kNone, kBlue, kOpenSquare, 1.0);
  SetHistProps(hMeasI, kBlue-1, kNone, kBlue-1, kDot, 1.0);
  hMeasI->SetLineWidth(2);

  UnfoldingUtils uu(hResp, hMeas, 0, 0, hTrue);

  // SVD analysis plot -----------------------------------------------
  // -----------------------------------------------------------------
  // Note that the Fourier coeffs. U'*b level off at deltab.
  SVDResult svdhists = uu.SVDAnalysis();
  uu.DrawSVDPlot(svdhists, 1e-5, 1e10);
  DrawObject(svdhists.U, "surf");

  // GSVD analysis ---------------------------------------------------
  // -----------------------------------------------------------------
  TMatrixD L = uu.LMatrix(n, UnfoldingUtils::k2DerivBC0);
  GSVDResult gsvd = uu.GSVDAnalysis(L, 0.38);
  DrawObject(gsvd.UHist, "surf");
  DrawObject(gsvd.UTbAbs,"pl");
  gsvd.UTbAbs->GetYaxis()->SetRangeUser(1e-5, 1e5);
  gPad->SetLogy();
  gsvd.coeffAbs->Draw("plsame");
  gsvd.regcAbs->Draw("plsame");

  DrawObject(hMeas, "pl", "Shaw test problem;x", cList);
  hMeasI->Draw("plsame");
  hTrue->Draw("plsame");
  gsvd.xregHist->Draw("plsame");

  // General-form Tikhonov algorithm using GSVD ----------------------
  // -----------------------------------------------------------------
  int nLambda = 80;
  TVectorD regVector(nLambda);
  for (int k=0; k<nLambda; k++) 
    regVector(k) = 0.02*(k+1);

  UnfoldingResult rg = uu.UnfoldTikhonovGSVD(gsvd, regVector);
  rg.XRegHist->SetFillColor(kCyan);
  rg.XRegHist->Draw("surf");

  DrawObject(rg.GcvCurve, "alp");
  SetGraphProps(rg.GcvCurve,kMagenta+2,kMagenta+2,kFullCircle,0.5);
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{min} = %g at k = %d", 
			      rg.lambdaGcv, rg.kGcv));
  TGraph* ggcv = new TGraph(1); 
  ggcv->SetPoint(0,rg.lambdaGcv,rg.GcvCurve->GetY()[rg.kGcv]);
  SetGraphProps(ggcv,kRed,kRed,kOpenCircle,2);
  ggcv->SetLineWidth(2);
  ggcv->Draw("psame");

  DrawObject(rg.LCurve, "alp");
  SetGraphProps(rg.LCurve,kBlue,kBlue,kFullCircle,0.5);

  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nIterPCGLS = 7;
  int LMatrixType = UnfoldingUtils::k2DerivBC0;
  uu.UnfoldPCGLS(nIterPCGLS,histsPCGLS,extrasPCGLS,LMatrixType);

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 500;
  TH1D* hX0 = hMeas; // Initial guess or "prior"
  uu.UnfoldRichardsonLucy(nIterRL, "", hX0);
  
  // SVD method (A. Hocker) ------------------------------------------
  // -----------------------------------------------------------------
  TString svdBC = "BC0, ~"; // Favor x=0 at edges (Reflect with "BCR")
  double lambda = 5.0;
  hSVD = uu.UnfoldSVD(lambda, genSvdAna, svdBC);
  if (0) { // Scan lambda regularization values
    double stepSize = 0.2;
    for (int i=0; i<50; i++) {
      lambda = stepSize*i;
      TH1D* h = uu.UnfoldSVD(lambda, svdResid, svdBC);
      svdAnim->Add(h);
    }
    TGraphTime* an = Animation(svdAnim, statObjs, "pl", 100 /*ms*/);
    DrawObject(an, "", "SVD", cList, 700, 500);
  }

  //  uu.DrawSVDPlot(svdHists, 1e-18, 1e18);
  uu.DrawGSVDPlot(genSvdAna, 1e-5, 1e2);

  if (0) {
    TGraph* svdRes = uu.ResidualNorm(svdResid, stepSize);
    double errNorm = svdRes->GetY()[0] + 2*TMath::Sqrt(n)*deltab;
    TLine e2Line;
    DrawObject(svdRes, "ALP");
    e2Line.DrawLine(0, errNorm, stepSize*svdRes->GetN(), errNorm);
    //  svdRes->GetYaxis()->SetRangeUser(0., 6);
  }

  // TH2D* hWcov = (TH2D*)genSvdAna->FindObject("hWcov1");
  // DrawObject(hWcov, "colz");
  // TH2D* hXcov = (TH2D*)genSvdAna->FindObject("hXcov1");
  // DrawObject(hXcov, "colz");


  //  return;

  DrawAll();
}

void DrawAll()
{
  DrawObject(hResp, "colz", "Response matrix;x_{meas};x_{true}", cList, 550, 500);
  gPad->SetRightMargin(0.15);
  
  // PCGLS hist properties
  for (int i=0; i<histsPCGLS->GetEntries(); i++) {
    TH1D* h = (TH1D*)histsPCGLS->At(i);
    SetHistProps(h, kRed, kNone, kRed, kOpenCircle, 1.5);
  }
  // Richardson-Lucy hist properties
  // for (int i=0; i<histsRL->GetEntries(); i++) {
  //   TH1D* h = (TH1D*)histsRL->At(i);
  //   SetHistProps(h, kRed+2, kNone, kRed+2, kOpenCircle, 1.5);
  // }
  
  // L-Curves  
  TGraph* pcgLCurve = (TGraph*)extrasPCGLS->At(0);
  DrawObject(pcgLCurve, "alp", "", cList);

  TLatex lt;
  lt.SetTextColor(kRed);
  for (int k=0; k<pcgLCurve->GetN(); k++) {
    double x = pcgLCurve->GetX()[k];
    double y = pcgLCurve->GetY()[k];
    lt.DrawLatex(x, y, Form("%d", k+1)); 
  }

  //  TGraph* rlLCurve = (TGraph*)extrasRL->At(0);
    // gL->SetLineColor(kRed);
    // gL->SetMarkerColor(kRed);
    // gL->SetMarkerStyle(kFullCircle);
    // gL->SetMarkerSize(1.5);
    // gL->SetLineWidth(2);

  //  DrawObject(rlLCurve, "alp", "", cList);

  TGraphTime* anim1 = Animation(histsPCGLS, statObjs, "pl", 200 /*ms*/);
  DrawObject(anim1, "1", "PCGLS", cList, 700, 500);

  // TGraphTime* anim2 = Animation(histsRL, statObjs, "pl", 0 /*ms*/);
  // DrawObject(anim2, "2", "Richardson-Lucy", cList, 700, 500);
 

  // Best solutions
  DrawObject(hMeas, "pl", "Shaw test problem;x", cList);
  hMeasI->Draw("plsame");
  hTrue->Draw("plsame");
  hSVD->SetLineWidth(2);
  hSVD->Draw("lsame");
  //  ((TH1D*)histsRL->At(499))->Draw("psame");
  ((TH1D*)histsPCGLS->At(4))->Draw("psame");
}

void LoadLibs()
{
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
}
