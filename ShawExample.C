bool savePDF = 0;
const int n=32;                     // number of points

TH1D *hTrue=0;
TH1D *hMeas=0, *hMeasI=0;       // Measured vector b, w/ & w/o noise
TH2D *hResp=0, *hRespT=0;       // Response matrix & transpose
TH1D* hSVD=0;

TObjArray* cList = new TObjArray();
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

  uu.DrawGSVDPlot(gsvd, 1e-5, 1e2);

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
  UnfoldingResult cg = uu.UnfoldPCGLS(nIterPCGLS,LMatrixType);
  DrawObject(cg.XRegHist,"surf");

  SetGraphProps(cg.LCurve,kRed,kRed,kFullCircle,1.2);
  DrawObject(cg.LCurve, "alp");
  TLatex ltx;
  ltx.SetTextColor(kRed);
  for (int k=0; k<cg.LCurve->GetN(); k++) {
    double x = cg.LCurve->GetX()[k];
    double y = cg.LCurve->GetY()[k];
    ltx.DrawLatex(x, y, Form("%d", k+1)); 
  }

  TGraphTime* anim1 = Animation(cg.XRegHist, statObjs, "pl", 200 /*ms*/, 
				kRed, kOpenCircle);
  DrawObject(anim1, "1", "PCGLS", cList, 700, 500);
  
  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 500;
  TH1D* hX0 = hMeas; // prior
  UnfoldingResult rl = uu.UnfoldRichardsonLucy(nIterRL, "", hX0);
  DrawObject(rl.XRegHist,"surf");
  SetGraphProps(rl.LCurve,kRed+2,kRed+2,kFullCircle,0.5);
  DrawObject(rl.LCurve, "alp", "", cList);
  TGraphTime* anim2 = Animation(rl.XRegHist, statObjs, "pl", 0,
				kRed+2,kFullCircle);
  DrawObject(anim2, "2", "Richardson-Lucy", cList, 700, 500);
  
  // SVD method (A. Hocker) ------------------------------------------
  // -----------------------------------------------------------------
  TString svdBC = "BC0,~"; // Favor x=0 at edges (Reflect with "BCR")
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

  DrawObject(hResp, "colz", "Response matrix;x_{meas};x_{true}", cList, 550, 500);
  gPad->SetRightMargin(0.15);

  // Best solutions
  DrawObject(hMeas, "pl", "Shaw test problem;x", cList);
  hMeasI->Draw("plsame");
  hTrue->Draw("plsame");
  hSVD->Draw("lsame");
  rg.hGcv->Draw("psame");
  cg.hLbest = cg.XRegHist->ProjectionX(Form("cg%d",4),4,4);
  cg.hLbest->Draw("psame");
  SetHistProps(rg.hGcv, kGreen+2, kNone, kGreen+2, kFullCircle, 1.0);
  SetHistProps(cg.hLbest, kRed, kNone, kRed, kOpenCircle, 1.0);

  return;
}

void LoadLibs()
{
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
}
