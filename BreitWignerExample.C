// Crosscheck with TSVDunfold and RooUnfold
// See $ROOTSYS/tutorials/math/TSVDUnfoldExample.C

Int_t nbins = 40;
TH1D *xini = new TH1D("xini", "MC truth", nbins, -10.0, 10.0);
TH1D *bini = new TH1D("bini", "MC reco", nbins, -10.0, 10.0);
TH2D *Adet = new TH2D("Adet", "detector response", nbins, -10.0, 10.0, nbins, -10.0, 10.0);
TH1D *data = new TH1D("data", "data", nbins, -10.0, 10.0);
TH1D *datatrue = new TH1D("datatrue", "data truth", nbins, -10.0, 10.0);
TH2D *statcov = new TH2D("statcov", "covariance matrix", nbins, -10.0, 10.0, nbins, -10.0, 10.0);
TRandom3 R;
const Double_t cutdummy= -99999.0;

// Unfolding solutions
TH1D* hSVD=0;

// Output containers
TObjArray* cList = new TObjArray();
TObjArray* histsRL = new TObjArray();
TObjArray* extrasRL = new TObjArray();
TObjArray* histsPCGLS = new TObjArray();
TObjArray* extrasPCGLS = new TObjArray();
TObjArray* svdHists = new TObjArray();
TObjArray* gsvdHists = new TObjArray();
TObjArray* statObjs = new TObjArray();

void BreitWignerExample()
{
  gStyle->SetOptStat(0);
  Load();

  SetHistProps(datatrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(data, kBlue, kNone, kBlue, kOpenSquare, 1.0);

  // Data/MC toy generation
  FillDistributions();
  UnfoldingUtils uu(Adet, data, 0, xini, datatrue, 0);

  // SVD analysis of the system
  uu.SVDAnalysis(svdHists);

  double lambda = 50.0;
  hSVD = uu.UnfoldSVD(lambda, gsvdHists, "BC0, ~");

  // Draw
  SetHistProps(hSVD, kGreen+2, kNone, kGreen+2, kOpenCircle, 1.0);

  DrawObject(Adet,"colz");
  DrawObject(datatrue, "ep", "Breit-Wigner test problem");
  data->Draw("epsame");
  hSVD->Draw("same");

  uu.DrawSVDPlot(svdHists, 1e-4, 1e4);
  uu.DrawGSVDPlot(gsvdHists, 5e-3, 100);

  // Covariance matrices
  TH2D* hWcov = (TH2D*)gsvdHists->FindObject("hWcov1");
  DrawObject(hWcov, "colz");
  TH2D* hXcov = (TH2D*)gsvdHists->FindObject("hXcov1");
  DrawObject(hXcov, "colz");
  TH2D* hXinv = (TH2D*)gsvdHists->FindObject("hXinv1");
  DrawObject(hXinv, "colz");
  TH2D* x1x = uu.TH2Product(hXinv, hXcov, "x1x");
  TH2D* xx1x = uu.TH2Product(hXcov, x1x, "XXinvX"); // supposed to equal hXcov
  DrawObject(xx1x, "colz");

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
