int m = 25, n = 25;
double x1 = 0., x2 = 5.;
double d = (x2-x1)/m;

// Response matrix, truth model, truth model estimator, and data.
TH2D* hResp=0;
TH1D* hTrue=0;
TH1D* hTrueData=0;
TH1D* hMeas=0;

// Unfolding solutions
TH1D* hSVD=0;
TH1D* hRL=0;
TH1D* hCh2=0;

TObjArray* cList = new TObjArray(); // Canvas container

void ConvolutionExample()
{
  gStyle->SetBarWidth(0.8);
  gStyle->SetOptStat(0);
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");

  // Create response matrix
  UnfoldingUtils utils;

  TF1* fTrue = new TF1("fTrue", "TMath::Exp(-x)", x1,x2);
  TF1* fGaus = new TF1("fGaus", "TMath::Gaus(x)", -10, 10);

  TestProblem t = 
    utils.MonteCarloConvolution(m, n, x1, x2, x1, x2, fTrue, fGaus, 10000);
  
  hResp = t.Response;
  hTrue = t.xTruth;
  hTrueData = t.xTruthEst;
  hMeas = t.bNoisy;

  // Create UnfoldingUtils instance
  // Use measured data as guess for solution. Clone possible only for m = n.
  TH2D* hMeasCov = 0;
  TH1D* hEff  = hResp->ProjectionY("hEff");
  UnfoldingUtils uu(hResp, hMeas, hMeasCov, 0, 0, hEff);

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 200;
  hRL = uu.UnfoldRichardsonLucy(nIterRL);

  // Chi squared minimization ----------------------------------------
  // -----------------------------------------------------------------
  uu.SetRegType(UnfoldingUtils::kTotCurv);
  hCh2 = uu.UnfoldChiSqMin(5e-4);

  // SVD algorithm ---------------------------------------------------
  // -----------------------------------------------------------------
  double lambda = 0.01;
  TObjArray* svdHists = new TObjArray();
  hSVD = uu.UnfoldSVD(lambda, svdHists, "BCR~");

  // -----------------------------------------------------------------
  // Draw
  // -----------------------------------------------------------------
  SetHistProps(hMeas, kBlue, kNone, kBlue, kOpenSquare, 0.8);
  SetHistProps(hTrue, kBlack, kNone, kBlack, kFullCircle, 0.8);
  SetHistProps(hTrueData, kBlack, kNone, kBlack, kFullCircle, 0.8);
  hTrue->SetLineWidth(2);
  SetHistProps(hSVD, kGreen+2, kNone, kGreen+2, kOpenCircle, 1.0);
  // SetHistProps(hCG, kRed, kNone, kRed, kOpenCircle, 1.5);
  SetHistProps(hCh2, kMagenta+2, kNone, kMagenta+2, kFullCircle, 1.0);
  SetHistProps(hRL, kRed+2, kNone, kRed+2, kOpenCircle, 1.2);

  // Do SVD analysis, put results in hists array
  TObjArray* hists = new TObjArray();
  uu.SVDAnalysis(hists, hResp, hMeas, "U");

  // Draw response matrix
  DrawObject(hResp, "col", "", cList, 500, 500);

  DrawObject(hEff, "", "Efficiency;t;#epsilon_{j}", cList, 700, 500);
  DrawObject(hSVD);

  // Draw true & meas
  DrawObject(hTrue, "l", "", cList, 700, 500);
  //  hTrueData->Draw("epsame");
  hMeas->Draw("epsame");
  //  hSVD->Draw("epsame");
  hRL->Draw("epsame");
  hCh2->Draw("epsame");
  TLegend* l1 = new TLegend(0.5, 0.6, 0.99, 0.99);
  l1->SetFillColor(kNone);
  l1->AddEntry(hTrue, "Data input model", "l");
  l1->AddEntry(hMeas, "Measurement", "epl");
  l1->AddEntry(hRL, Form("Richardson-Lucy alg."), "epl");
  l1->AddEntry(hCh2, Form("#chi^{2} minimization method", nIterRL), "epl");
  l1->Draw();

  // Draw a few of the left singular vectors
  TH1D* hu[999];
  for (int i=0; i<n; i++) {
    hu[i] = (TH1D*)hists->FindObject(Form("SV_u_%d",i));
    hu[i]->SetLineWidth(2);
    hu[i]->SetTitle(Form("Left singular vector u_{%d};Column index",i));
  }
  SetHistProps(hu[n-1], kGray, kNone, kGray, kFullCircle, 0.8);
  DrawObject(hu[n-1], "pl", 
	     Form("Left Singular Vectors "
		  "#color[4]{u_{1}}, "
		  "#color[2]{u_{2}}, "
		  "#color[418]{u_{3}}, "
		  "#color[922]{u_{%d}}", n),
  	     cList);
  SetHistProps(hu[2], kGreen+1, kNone, kGreen+1, kFullCircle, 0.8);
  hu[2]->Draw("plsame");
  SetHistProps(hu[1], kRed, kNone, kRed, kFullCircle, 0.8);
  hu[1]->Draw("plsame");
  SetHistProps(hu[0], kAzure-2, kNone, kAzure-2, kFullCircle, 0.8);
  hu[0]->Draw("plsame");

  // Draw the s.v. spectrum
  TH1D* hsig = hists->FindObject("hsig1");
  hsig->SetTitle("Singular values;Column index;#sigma_{i}");
  SetHistProps(hsig, kBlack, kNone, kBlack, kFullCircle, 0.8);
  hsig->SetLineWidth(2);
  DrawObject(hsig, "p", "", cList, 500, 500);
  gPad->SetLogy();

  // Draw s.v. coefficients |u_i'*b|
  cList->Add(uu.DrawSVDPlot(hists, 0.2, 1e9));

  if (1)
    PrintPDFs(cList, "./outputs");

}
