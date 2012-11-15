bool savePDF = 0;
const int n=32;                     // number of points

TH1D *hTrue=0, *hXini=0;
TH1D *hMeas=0, *hMeasI=0;       // Measured vector b, w/ & w/o noise
TH2D *hResp=0, *hRespT=0;       // Response matrix & transpose
TH1D* hSVD=0;

TObjArray* cList = new TObjArray();
TObjArray* histsRL = new TObjArray();
TObjArray* extrasRL = new TObjArray();
TObjArray* histsPCGLS = new TObjArray();
TObjArray* extrasPCGLS = new TObjArray();
TObjArray* svdHists = new TObjArray();
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

  // Setup:
  UnfoldingUtils shaw;
  shaw.SetTrueRange(0.,1); 
  shaw.SetMeasRange(0.,1);
  double deltab = 0.05; // Gaussian white noise on b (absolute)

  TestProblem p = shaw.ShawSystem(n,deltab);
  hResp  = p.Response;
  hTrue  = p.xTruth;
  hMeas  = p.bNoisy;
  hMeasI = p.bIdeal;
  hXini  = hResp->ProjectionX("hXini");

  statObjs->Add(hMeasI); 
  statObjs->Add(hMeas); 
  statObjs->Add(hTrue);
  SetHistProps(hTrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(hMeas, kBlue, kNone, kBlue, kOpenSquare, 1.0);
  SetHistProps(hMeasI, kBlue-1, kNone, kBlue-1, kDot, 1.0);
  hMeasI->SetLineWidth(2);

  //  hXini->Scale(hTrue->Integral()/hXini->Integral());
  for (int i=0; i<n; i++)
    hXini->SetBinContent(i+1,1.0);


  UnfoldingUtils uu(hResp, hMeas, 0, hXini, hTrue);

  // SVD Picard plot -------------------------------------------------
  // -----------------------------------------------------------------
  uu.SVDAnalysis(svdHists);

  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nIterPCGLS = 7;
  uu.UnfoldPCGLS(nIterPCGLS,histsPCGLS,extrasPCGLS,UnfoldingUtils::k2DerivBC0, "~");

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 500;
  TH1D* hX0 = hMeas; // Initial guess or "prior"
  uu.UnfoldRichardsonLucy(nIterRL, histsRL, extrasRL, "", hX0);
  
  // SVD method (A. Hocker) ------------------------------------------
  // -----------------------------------------------------------------
  TString svdBC = "BC0, ~"; // Favor x=0 at edges (Reflect with "BCR")
  double lambda = 5.0;
  hSVD = uu.UnfoldSVD(lambda, genSvdAna, svdBC);
  if (1) { // Scan lambda regularization values
    double stepSize = 0.1;
    for (int i=0; i<100; i++) {
      lambda = stepSize*i;
      TH1D* h = uu.UnfoldSVD(lambda, svdResid, svdBC);
      svdAnim->Add(h);
    }
    TGraphTime* an = Animation(svdAnim, statObjs, "pl", 100 /*ms*/);
    DrawObject(an, "", "SVD", cList, 700, 500);
  }

  uu.DrawSVDPlot(svdHists, 1e-18, 1e18);
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
  for (int i=0; i<histsRL->GetEntries(); i++) {
    TH1D* h = (TH1D*)histsRL->At(i);
    SetHistProps(h, kRed+2, kNone, kRed+2, kOpenCircle, 1.5);
  }
  
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

  TGraph* rlLCurve = (TGraph*)extrasRL->At(0);
  DrawObject(rlLCurve, "alp", "", cList);

  TGraphTime* anim1 = Animation(histsPCGLS, statObjs, "pl", 500 /*ms*/);
  DrawObject(anim1, "1", "PCGLS", cList, 700, 500);

  TGraphTime* anim2 = Animation(histsRL, statObjs, "pl", 0 /*ms*/);
  DrawObject(anim2, "2", "Richardson-Lucy", cList, 700, 500);
 

  // Best solutions
  DrawObject(hMeas, "pl", "Shaw test problem;x", cList);
  hMeasI->Draw("plsame");
  hTrue->Draw("plsame");
  hSVD->SetLineWidth(2);
  hSVD->Draw("lsame");
  ((TH1D*)histsRL->At(499))->Draw("psame");
  ((TH1D*)histsPCGLS->At(5))->Draw("psame");
}

void LoadLibs()
{
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
}
