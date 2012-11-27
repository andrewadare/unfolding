// --- Apply UnfoldingUtils to RooUnfold test problem ---
//
// See $ROOTSYS/tutorials/math/TSVDUnfoldExample.C or
// RooUnfoldExample.cxx in the RooUnfold package.
//
// The MC response matrix is generated from a Breit-Wigner truth
// distribution. The "data" is Gaussian. Both are distorted by
// Reconstruct(), which includes the following:
// - systematic offset bias (2.5 units)
// - variable efficiency (linear, 30-100%)
// - modest Gaussian smearing (sigma = 0.2).
//
// Four unfolding methods are applied from UnfoldingUtils.

Int_t nbins = 40;           // square problem (m = n)
double xt1 = -10., xt2=10.; // range of true variable
double xm1 = -10., xm2=10.; // range of meas variable

TH1D *xini = new TH1D("xini", 
		      "Black: MC truth (Breit-Wigner), "
		      "Blue: reconstructed;x;counts",nbins,xt1,xt2);
TH1D *bini = new TH1D("bini","MC reco",nbins,xm1,xm2);
TH2D *Adet = new TH2D("Adet","detector response",
		      nbins,xm1,xm2,nbins,xt1,xt2);
TH1D *data = new TH1D("data","data",nbins,xm1,xm2);
TH1D *datatrue = new TH1D("datatrue", 
			  "Black: data truth (Gaussian), "
			  "Blue: reconstructed;x;counts",
			  nbins,xm1,xm2);
TH2D *statcov = new TH2D("statcov", "covariance matrix", 
			 nbins,xm1,xm2,nbins,xt1,xt2);
TH1D* hxeff = (TH1D*)xini->Clone("hxeff");

TLatex lt;
TRandom3 R;
const Double_t cutdummy= -99999.0;

// Unfolding solutions
TH1D* hSVD=0;
TH1D* hCh2=0;
TH1D* hCG=0;
TH1D* hRL=0;

// Output containers
TObjArray* cList = new TObjArray();
// TObjArray* histsRL = new TObjArray();
// TObjArray* extrasRL = new TObjArray();
// TObjArray* histsCG = new TObjArray();
// TObjArray* extrasCG = new TObjArray();
// TObjArray* svdHists = new TObjArray();
TObjArray* gsvdHists = new TObjArray();
TObjArray* statObjs = new TObjArray();

void RUExample()
{
  gStyle->SetOptStat(0);
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
  
  // Fill hxeff with efficiency from Reconstruct()
  hxeff->SetTitle("Efficiency;x;efficiency");
  for (int j=1; j<=nbins; j++) {
    double xt = hxeff->GetBinCenter(j);
    Double_t eff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);
    hxeff->SetBinContent(j,eff);
  }

  // Data/MC toy generation
  FillDistributions();
  UnfoldingUtils uu(Adet, data, 0, xini, datatrue, hxeff);

  // SVD and GSVD analysis of the system -----------------------------
  // -----------------------------------------------------------------
  SVDResult svd = uu.SVDAnalysis(0,0,"~");
  uu.DrawSVDPlot(svd, 1e-4, 1e4);

  TMatrixD L = uu.LMatrix(nbins, UnfoldingUtils::k2DerivBC0);
  GSVDResult gsvd = uu.GSVDAnalysis(L,100,0,0,"~");
  uu.DrawGSVDPlot(gsvd, 1e-4, 2e2);

  // General-form Tikhonov algorithm using GSVD ----------------------
  // -----------------------------------------------------------------
  int nLambda = 50;
  TVectorD regVector(nLambda);
  for (int k=0; k<nLambda; k++) 
    regVector(k) = k+50;
  
  UnfoldingResult rg = uu.UnfoldTikhonovGSVD(gsvd, regVector);
  DrawObject(rg.XRegHist, "surf");
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

  // SVD algorithm ---------------------------------------------------
  // -----------------------------------------------------------------
  double lambda = 50.0;
  hSVD = uu.UnfoldSVD(lambda, gsvdHists, "BC0, ~");

  // Chi squared minimization ----------------------------------------
  // -----------------------------------------------------------------
  TVectorD regWts(40);
  for (int k=0; k<40; k++)
    regWts(k) = (k+1)*2e-5;
  uu.SetRegType(UnfoldingUtils::kTotCurv);
  UnfoldingResult cs = uu.UnfoldChiSqMin(regWts);
  DrawObject(cs.LCurve,"alp");
  SetGraphProps(cs.LCurve,kBlue,kBlue,kFullCircle,0.5);
  hCh2 = cs.XRegHist->ProjectionX(Form("cs%d",20),20,20);

  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nCG = 9;
  UnfoldingResult cg = 
    uu.UnfoldPCGLS(nCG, UnfoldingUtils::k2DerivBC0,"~");
  DrawObject(cg.LCurve,"alp");
  SetGraphProps(cg.LCurve,kBlue,kBlue,kFullCircle,0.5);
  hCG = cg.XRegHist->ProjectionX(Form("cg%d",5),5,5);
  
  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 10;
  TH1D* hX0 = data->Clone("hX0"); // use measurement as prior
  hX0->Smooth(5);
  UnfoldingResult rl = uu.UnfoldRichardsonLucy(nIterRL, "^", hX0);
  hRL = rl.XRegHist->ProjectionX(Form("rl%d",nIterRL),nIterRL,nIterRL);

  // -----------------------------------------------------------------
  // Draw
  // -----------------------------------------------------------------
  SetHistProps(datatrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(data, kBlue, kNone, kBlue, kOpenSquare, 1.0);
  SetHistProps(hSVD, kGreen, kNone, kGreen, kOpenCircle, 1.0);
  SetHistProps(rg.hGcv, kGreen+2, kNone, kGreen+2, kOpenCircle, 1.0);
  SetHistProps(hCG, kRed, kNone, kRed, kOpenCircle, 1.5);
  SetHistProps(hCh2, kMagenta+2, kNone, kMagenta+2, kFullCircle, 1.0);
  SetHistProps(hRL, kRed+2, kNone, kRed+2, kOpenCircle, 1.2);

  // MC and data histos
  DrawObject(hxeff);
  DrawObject(xini);
  bini->SetLineColor(kBlue);
  bini->Draw("same");
  DrawObject(datatrue);
  data->Draw("same");

  // Response matrix, as generated (counts)
  DrawObject(Adet,"colz");

  // Probability resp. matrix (normalized to efficiency)
  TH2D* hAprob = uu.GetAProbHist();
  hAprob->SetTitle("Detector response (normalized to efficiency);"
		   "measured;true");
  DrawObject(hAprob,"colz");

  // Unfolding results
  DrawObject(datatrue, "ep");
  data->Draw("epsame");
  hSVD->Draw("epsame");
  rg.hGcv->Draw("epsame");
  hCh2->Draw("epsame");
  hCG->Draw("epsame");
  hRL->Draw("epsame");

  /*
  // Covariance matrices
  if (0) {
    // Regularization covariance from SVD method
    TH2D* hXcov = (TH2D*)gsvdHists->FindObject("hXcov1");
    DrawObject(hXcov, "colz");

    // Statistical error propagation via 100 MC pseudo experiments.
    // Works for any unfolding method, SVD shown here.
    TH2D* hUcov = uu.UnfoldCovMatrix(100, UnfoldingUtils::kSVDAlgo, 
				     lambda, "BC0,~");
    DrawObject(hUcov, "colz");
    hXcov->Add(hUcov); // Add independent contributions
    for (int j=0; j<hSVD->GetNbinsX(); j++)
      hSVD->SetBinError(j, TMath::Sqrt(hXcov->GetBinContent(j,j)));
  }
  */

  // Animations
  statObjs->Add(datatrue);
  statObjs->Add(data);

  TGraphTime* aRL = Animation(rl.XRegHist, statObjs, "ep", 100, kRed+2);
  TGraphTime* aCG = Animation(cg.XRegHist, statObjs, "ep", 100, kRed);

  if (0) {
    DrawObject(aRL, "2", "Richardson-Lucy", cList, 700, 500);
    DrawObject(aCG, "2", "Conjugate Gradient for Least Squares", cList, 700, 500);
  }
  else {
    // Avoid appending to old gifs
    TString cmd("for i in *.gif; do rm $i; done");
    gSystem->Exec(cmd.Data());

    TCanvas* crl = new TCanvas("crl", "crl", 1);
    aRL->SaveAnimatedGif("RichLucyAnim.gif");

    TCanvas* ccg = new TCanvas("ccg", "ccg", 1);
    aCG->SaveAnimatedGif("PCGLSAnim.gif");

    // Convert gif --> pdf
    cmd = "for i in *.gif; do convert $i ${i%.*}.pdf; done";
    gSystem->Exec(cmd.Data());
  }

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
