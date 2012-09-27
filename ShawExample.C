bool savePDF = 0;
const int n=32;                     // number of points

TH1D *hTrue=0, *hXini=0;
TH1D *hMeas=0, *hMeasIdeal=0;       // Measured vector b w/0, w/ noise
TH2D *hResp=0, *hRespT=0;           // Response matrix & transpose
TH1D* hSVD=0;

TObjArray* cList = new TObjArray();
TObjArray* histsRL = new TObjArray();
TObjArray* extrasRL = new TObjArray();
TObjArray* histsPCGLS = new TObjArray();
TObjArray* extrasPCGLS = new TObjArray();
TObjArray* svdHists = new TObjArray();
TObjArray* gsvdHists = new TObjArray();
TObjArray* statObjs = new TObjArray();

TLatex lt;

void ShawExample()
{
  gStyle->SetOptStat(0);
  lt.SetNDC();
  LoadLibs();
  UnfoldingUtils uu;
  uu.SetTrueRange(0.,1); 
  uu.SetMeasRange(0.,1);
  double deltab = 0.05; // Gaussian white noise on b

  // Setup:
  TObjArray* ShawArray = uu.ShawSystem(n,deltab);
  hResp = (TH2D*)ShawArray->FindObject("Shaw_A");
  hTrue = (TH1D*)ShawArray->FindObject("Shaw_x");
  hMeas = (TH1D*)ShawArray->FindObject("Shaw_b");
  hXini = hResp->ProjectionX("hXini");
  statObjs->Add(hMeas); 
  statObjs->Add(hTrue);
  SetHistProps(hTrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(hMeas, kBlue, kNone, kBlue, kOpenSquare, 1.0);

  // SVD Picard plot -------------------------------------------------
  // -----------------------------------------------------------------
  uu.SVDAnalysis(hResp, hMeas, svdHists);

  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nIterPCGLS = 7;
  uu.UnfoldPCGLS(hResp,hMeas,0,nIterPCGLS,histsPCGLS,extrasPCGLS,
		 UnfoldingUtils::k2DerivBC0);

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 500;
  TH1D* hX0 = hXini; // Initial guess or "prior"
  uu.UnfoldRichardsonLucy(hResp, hMeas, hX0, nIterRL, histsRL, 
			  extrasRL, hXini);
  
  // SVD method (A. Hocker) ------------------------------------------
  // -----------------------------------------------------------------
  TString svdBC = "BC0"; // Favor x=0 at edges (Reflect with "BCR")
  double lambda = 50;
  hSVD = uu.UnfoldSVD(hResp, hMeas, hXini, lambda, gsvdHists, svdBC);
  if (1) { // Scan lambda regularization values
    for (int i=0; i<20; i++) {
      lambda = 5*i;
      TH1D* h = uu.UnfoldSVD(hResp, hMeas, hXini, lambda, 0, svdBC);
      gsvdHists->Add(h);
    }
    TGraphTime* an = Animation(gsvdHists, statObjs, "pl", 100 /*ms*/);
    DrawObject(an, "", "SVD", cList, 700, 500);
  }

  uu.DrawSVDPlot(svdHists, 1e-18, 1e18);
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
