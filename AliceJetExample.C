bool savePDF = 0;
TDatime dtime;
TLatex lt;
TString outFileName = Form("outputs/PbPbSpectra-Unfolded-%d.root", dtime.GetDate());
TString pdfName = Form("plots/PbPb-%d", dtime.GetDate()); // .pdf extension gets appended

TFile* inFile = new TFile("inputs/andrew-20120920.root", "read");
TFile* outFile = new TFile(outFileName, "recreate");
TH2D* hResp = (TH2D*)inFile->Get("hResp");
TH1D* hMeas = (TH1D*)inFile->Get("hMeas");
TH1D* hXini = (TH1D*)inFile->Get("hXini");

TH1D* hSVD = 0, *hChi2 = 0;

TObjArray* cList = new TObjArray();
TObjArray* svdHists = new TObjArray();
TObjArray* gsvdHists = new TObjArray();
TObjArray* histsPCGLS = new TObjArray();
TObjArray* extrasPCGLS = new TObjArray();
TObjArray* histsRL = new TObjArray();
TObjArray* extrasRL = new TObjArray();
TObjArray* statObjs = new TObjArray();

void UnfoldingMacro()
{
  gStyle->SetOptStat(0);
  lt.SetNDC();
  LoadLibs();

  SetHistProps(hMeas, kBlack, 0, kBlack, kFullCircle, 1.0);
  statObjs->Add(hMeas); 

  hXini->Scale(hMeas->Integral()/hXini->Integral());

  UnfoldingUtils uu(hResp, hMeas, hXini);
  Double_t xt1=uu.GetTrueX1(),xt2=uu.GetTrueX2(),
           xm1=uu.GetMeasX1(),xm2=uu.GetMeasX2();
  Int_t n = uu.GetN();               // true/unfolded bins
  Int_t m = uu.GetM();               // measured bins

  TH2D* hATilde = uu.GetATildeHist();
  TH1D* hbTilde = uu.GetbTildeHist();

  if (0) {
    // ================== TSVDUnfold (root) ============================
    // =================================================================
    TH1D* hbini = hResp->ProjectionX("hbini");
    TSVDUnfold *tsvdunf = new TSVDUnfold(hMeas, hbini, hXini, hResp);
    TH1D* unfresult = tsvdunf->Unfold(10);
    TH1D* sv = tsvdunf->GetSV();
    TH1D* d  = tsvdunf->GetD();
    SetHistProps(sv, kGreen+2, 0, kGreen+2, kOpenCircle, 1.0);
    SetHistProps(d,  kGreen+2, 0, kGreen+2, kFullSquare, 1.0);
    DrawObject(hMeas);
    hXini->Draw("same");
    unfresult->Draw("same");
    gPad->SetLogy();
    // return;
  }
  
  // ==================== Reg. Chi2 Min. =============================
  // =================================================================
  uu.SetRegType(UnfoldingUtils::kTotCurv);
  double regwt = 100.;
  hChi2 = uu.UnfoldChiSqMin(hResp, hMeas, 0, 0, hXini, regwt, 0);
  SetHistProps(hChi2, kGreen, kNone, kGreen, kFullSquare, 1.0);

  // ====================== (P)CGLS ==================================
  // =================================================================
  // Smoothing matrix L (dimension n x p) and its nullspace W
  Int_t nIterPCGLS = 12;
  uu.UnfoldPCGLS(hResp, hMeas, hXini, nIterPCGLS, histsPCGLS,   
		 extrasPCGLS, UnfoldingUtils::k2DerivBCR, "SB");
  
  // ==================== Richardson-Lucy ============================
  // =================================================================
  Int_t nIterRL = 50;
  // Use hXini dist. as the initial solution
  TH1D* hXStart = (TH1D*)hXini->Clone("hXStart");

  // Solve un-scaled problem
  uu.UnfoldRichardsonLucy(hResp, hMeas, hXStart, nIterRL, histsRL, extrasRL, 0);

  // =============== SVD (Hocker/Kartvelishvili) =====================
  // =================================================================
  // NIM A A372, 469 (1996) eq. (39)

  // SVD analysis of unregularized system
  uu.SVDAnalysis(hResp, hMeas, svdHists);

  double lambda = 100;
  hSVD = (TH1D*)uu.UnfoldSVD(hResp, hMeas, hXini, lambda, gsvdHists);

  uu.DrawSVDPlot(svdHists, 1e-5, 1e9); // full scale
  uu.DrawSVDPlot(svdHists, 1e3, 1e7); // zoom on SVD coeffs.
  uu.DrawGSVDPlot(gsvdHists, 1e-9, 1e5);

  // //   TSVD results
  // sv->Draw("plsame");
  // d->Draw("plsame");

  DrawAll();
  // TH2D* hZcov = (TH2D*)gsvdHists->FindObject("hZcov1");
  // DrawObject(hZcov, "colz");
  TH2D* hWcov = (TH2D*)gsvdHists->FindObject("hWcov1");
  DrawObject(hWcov, "colz");
  TH2D* hXcov = (TH2D*)gsvdHists->FindObject("hXcov1");
  DrawObject(hXcov, "colz");
  TH2D* hXinv = (TH2D*)gsvdHists->FindObject("hXinv1");
  DrawObject(hXinv, "colz");

  TH2D* x1x = uu.TH2Product(hXinv, hXcov, "x1x");
  TH2D* xx1x = uu.TH2Product(hXcov, x1x, "XXinvX");
  DrawObject(xx1x, "colz");
}

void DrawAll()
{
  // TH1D* respEff = hResp->ProjectionY("respEff");
  // DrawObject(respEff);

  DrawLCurves();
  DrawMeasAndUnfolded();

  // Richardson-Lucy animation
  TGraphTime* animRL = Animation(histsRL, statObjs, "pl", 10 /*ms*/);
  DrawObject(hMeas, "", "", cList, 700, 500);
  gPad->SetLogy();
  animRL->Draw("same");
  ((TH1D*)histsRL->At(4))->Draw("psame");
}

void DrawLCurves()
{
  TGraph* pcgLCurve   = (TGraph*)extrasPCGLS->At(0);
  SetGraphProps(pcgLCurve,   kRed,   kRed,   kFullCircle, 1.0);
  DrawObject(pcgLCurve, "alp", "", cList);
  TLatex lt;
  lt.SetTextSize(0.04);
  lt.SetTextColor(kRed);
  for (int k=0; k<pcgLCurve->GetN(); k++) {
    double x = pcgLCurve->GetX()[k];
    double y = pcgLCurve->GetY()[k];
    lt.DrawLatex(x, y, Form("%d", k+1)); 
  }

  TGraph* rlLCurve = (TGraph*)extrasRL->At(0);
  DrawObject(rlLCurve, "alp", "", cList);

  return;
}

void DrawMeasAndUnfolded()
{
  // ================ Measured and unfolded result =======================
  DrawObject(hMeas, "ep", "Measured and unfolded p_{T} distributions;"
	     "p_{T} (GeV/c);counts", cList, 700, 700);
  TAxis* ax = hMeas->GetXaxis();
  TAxis* ay = hMeas->GetYaxis();
  ax->CenterTitle();
  ax->SetTitleOffset(1.2);
  ay->SetTitleOffset(1.5);
  gPad->SetLogy();

  int nHistsCG_ = const_cast<int>(histsPCGLS->GetEntries());
  const int nHistsCG = nHistsCG_;
  TH1D* h[nHistsCG];
  for (int i=0; i<nHistsCG; i++) {
    h[i] = (TH1D*)histsPCGLS->At(i);
    SetHistProps(h[i], kBlue, kNone, kBlue, kOpenCircle, 1.0);
  }
  h[3]->Draw("epsame");

  // Richardson-Lucy
  int nHistsRL_ = const_cast<int>(histsRL->GetEntries());
  const int nHistsRL = nHistsRL_;
  TH1D* hrl[nHistsRL];
  for (int i=0; i<nHistsRL; i++) {
    hrl[i] = (TH1D*)histsRL->At(i);
    SetHistProps(hrl[i], kRed+2, kNone, kRed+2, kOpenCircle, 0.8);
  }
  hrl[4]->Draw("epsame");

  SetHistProps(hSVD, kGreen+2, kNone, kGreen+2, kFullCircle, 0.8);
  hSVD->Draw("epsame");

  hChi2->Draw("epsame");

  hXini->Draw("same");

  TLegend* leg = new TLegend(0.7, 0.65, 0.99, 0.9);
  leg->SetFillColor(kNone);
  leg->AddEntry(hMeas, "Data", "epl");
  leg->AddEntry(hXini, "x^{ini} (MC)", "epl");
  leg->AddEntry(h[3], "PCGLS", "epl");
  leg->AddEntry(hrl[0], "R-L", "epl");
  leg->AddEntry(hSVD, "SVD", "epl");
  leg->Draw();

  return;
}

void LoadLibs()
{
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UnfoldingUtils.C+g");
  gROOT->LoadMacro("IOUtilFns.C+g");
}
