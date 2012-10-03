// Crosscheck with TSVDunfold and RooUnfold
// See $ROOTSYS/tutorials/math/TSVDUnfoldExample.C

Int_t nbins = 40;
TH1D *xini = new TH1D("xini", "MC truth (Breit-Wigner)", nbins, -10.0, 10.0);
TH1D *bini = new TH1D("bini", "MC reco", nbins, -10.0, 10.0);
TH2D *Adet = new TH2D("Adet", "detector response", nbins, -10.0, 10.0, nbins, -10.0, 10.0);
TH1D *data = new TH1D("data", "data", nbins, -10.0, 10.0);
TH1D *datatrue = new TH1D("datatrue", "data truth", nbins, -10.0, 10.0);
TH2D *statcov = new TH2D("statcov", "covariance matrix", nbins, -10.0, 10.0, nbins, -10.0, 10.0);
TH1D* hxeff = (TH1D*)xini->Clone("hxeff");

TRandom3 R;
const Double_t cutdummy= -99999.0;

// Unfolding solutions
TH1D* hSVD=0;
TH1D* hCh2=0;
TH1D* hCG=0;
TH1D* hRL=0;

// Output containers
TObjArray* cList = new TObjArray();
TObjArray* histsRL = new TObjArray();
TObjArray* extrasRL = new TObjArray();
TObjArray* histsCG = new TObjArray();
TObjArray* extrasCG = new TObjArray();
TObjArray* svdHists = new TObjArray();
TObjArray* gsvdHists = new TObjArray();
TObjArray* statObjs = new TObjArray();

void BreitWignerExample()
{
  gStyle->SetOptStat(0);
  Load();

  for (int j=1; j<=nbins; j++) {
    double xt = hxeff->GetBinCenter(j);
    Double_t eff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);  // efficiency
    hxeff->SetBinContent(j,eff);
  }
  DrawObject(hxeff);

  // Data/MC toy generation
  FillDistributions();
  UnfoldingUtils uu(Adet, data, 0, xini, datatrue, hxeff);

  DrawObject(xini);
  bini->SetLineColor(kBlue);
  bini->Draw("same");

  DrawObject(datatrue);
  data->Draw("same");


  TH2D* hAprob = uu.GetAProbHist();
  hAprob->SetTitle("Detector response (normalized to efficiency);measured;true");
  DrawObject(hAprob,"colz");

  // SVD analysis of the system
  uu.SVDAnalysis(svdHists);

  // SVD algorithm ---------------------------------------------------
  // -----------------------------------------------------------------
  double lambda = 50.0;
  hSVD = uu.UnfoldSVD(lambda, gsvdHists, "BC0, ~");

  if (0) {
    // Covariance matrices
    TH2D* hWcov = (TH2D*)gsvdHists->FindObject("hWcov1");
    DrawObject(hWcov, "colz");
    TH2D* hXcov = (TH2D*)gsvdHists->FindObject("hXcov1");
    DrawObject(hXcov, "colz");
    TH2D* hXinv = (TH2D*)gsvdHists->FindObject("hXinv1");
    DrawObject(hXinv, "colz");
    TH2D* hUcov = uu.UnfoldCovMatrix(100, UnfoldingUtils::kSVDAlgo, lambda, "BC0,~");
    DrawObject(hUcov, "colz");
    hXcov->Add(hUcov); // Add independent contributions
    for (int j=0; j<hSVD->GetNbinsX(); j++)
      hSVD->SetBinError(j, TMath::Sqrt(hXcov->GetBinContent(j,j)));
  }

  // Chi squared minimization ----------------------------------------
  // -----------------------------------------------------------------
  uu.SetRegType(UnfoldingUtils::kTotCurv);
  double regwt = 5e-4;
  hCh2 = uu.UnfoldChiSqMin(regwt, 0, "");

  // regwt = 50;
  // hCh2 = uu.UnfoldChiSqMin(regwt, 0, "~");


  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nCG = 9;
  hCG = uu.UnfoldPCGLS(nCG, histsCG, extrasCG, 
		       UnfoldingUtils::k2DerivBC0, "~",0,0,xini);

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
    // for (int j=0; j<hSVD->GetNbinsX(); j++)
    //   hSVD->SetBinError(j, TMath::Sqrt(hXcov->GetBinContent(j,j)));

  int nIterRL = 8;
  TH1D* hX0 = data->Clone("hX0"); // Initial guess or "prior"
  hX0->Smooth(5);
  hRL = uu.UnfoldRichardsonLucy(nIterRL, histsRL, extrasRL, "^", hX0);

  // Draw
  SetHistProps(datatrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(data, kBlue, kNone, kBlue, kOpenSquare, 1.0);
  SetHistProps(hSVD, kGreen+2, kNone, kGreen+2, kOpenCircle, 1.0);
  SetHistProps(hCG, kRed, kNone, kRed, kOpenCircle, 1.5);
  SetHistProps(hCh2, kMagenta+2, kNone, kMagenta+2, kFullCircle, 1.0);
  SetHistProps(hRL, kRed+2, kNone, kRed+2, kOpenCircle, 1.2);

  // DrawObject(hCh2);
  // DrawObject(hRL);

  DrawObject(Adet,"colz");
  DrawObject(datatrue, "ep", "Data truth (Gaussian)");
  data->Draw("epsame");
  hSVD->Draw("epsame");
  hCh2->Draw("epsame");
  hCG->Draw("epsame");
  hRL->Draw("epsame");

  uu.DrawSVDPlot(svdHists, 1e-4, 1e4);
  uu.DrawGSVDPlot(gsvdHists, 5e-3, 100);

  statObjs->Add(datatrue);
  statObjs->Add(data);
  TGraphTime* aRL = Animation(histsRL, statObjs, "ep", 500 /*ms*/);
  DrawObject(aRL, "2", "Richardson-Lucy", cList, 700, 500);

  TGraphTime* aCG = Animation(histsCG, statObjs, "ep", 500 /*ms*/);
  DrawObject(aCG, "2", "Conjugate Gradient for Least Squares", cList, 700, 500);

}

void Load()
{
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
}

void FillDistributions()
{
  // Fill the MC using a Breit-Wigner, mean 0.3 and width 2.5.
  for (Int_t i= 0; i<100000; i++) {
    Double_t xt = R.BreitWigner(0.3, 2.5);
    xini->Fill(xt);
    Double_t x = Reconstruct( xt, R );
    if (x != cutdummy) {
      Adet->Fill(x, xt);
      bini->Fill(x);
    }
  }
  
  // Fill the "data" with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<10000; i++) {
    Double_t xt = R.Gaus(0.0, 2.0);
    datatrue->Fill(xt);
    Double_t x = Reconstruct( xt, R );
    if (x != cutdummy) 
      data->Fill(x);
  }
  
  // Fill the data covariance matrix
  for (int i=1; i<=data->GetNbinsX(); i++) {
    statcov->SetBinContent(i,i,data->GetBinError(i)*data->GetBinError(i)); 
  }
}

Double_t Reconstruct( Double_t xt, TRandom3& R )
{
  // apply some Gaussian smearing + bias and efficiency corrections to fake reconstruction
  const Double_t cutdummy = -99999.0;
  Double_t xeff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);  // efficiency
  Double_t x    = R.Rndm();
  if (x > xeff) return cutdummy;
  else {
    Double_t xsmear= R.Gaus(-2.5,0.2); // bias and smear
    return xt+xsmear;
  }
}
