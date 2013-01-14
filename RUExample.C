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

Int_t nb = 80;
Int_t nx = 40;
double xt1 = -10., xt2=10.; // range of true variable
double xm1 = -10., xm2=10.; // range of meas variable

TH1D *xini = new TH1D("xini", 
		      "Black: MC truth (Breit-Wigner), "
		      "Blue: reconstructed;x;counts",nx,xt1,xt2);
TH1D *bini = new TH1D("bini","MC reco",nb,xm1,xm2);
TH2D *Adet = new TH2D("Adet","detector response",
		      nb,xm1,xm2,nx,xt1,xt2);
TH1D *data = new TH1D("data","data",nb,xm1,xm2);
TH1D *datatrue = new TH1D("datatrue", 
			  "Black: data truth (Gaussian), "
			  "Blue: reconstructed;x;counts",
			  nx,xm1,xm2);
TH2D *statcov = new TH2D("statcov", "covariance matrix", 
			 nb,xm1,xm2,nx,xt1,xt2);
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
TObjArray* gsvdHists = new TObjArray();
TObjArray* statObjs = new TObjArray();

void RUExample()
{
  gStyle->SetOptStat(0);
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("UtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
  lt.SetNDC();

  // Fill hxeff with efficiency from Reconstruct()
  hxeff->SetTitle("Efficiency;x;efficiency");
  for (int j=1; j<=nx; j++) {
    double xt = hxeff->GetBinCenter(j);
    Double_t eff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);
    hxeff->SetBinContent(j,eff);
  }

  // Data/MC toy generation
  FillDistributions();
  TH1D* Aproj = Adet->ProjectionY("Aproj", 1,nb);

  Printf("Aproj: %g, xini %g, data %g, datatrue %g, eff %g",
	 Aproj->Integral(1,nx),
	 xini->Integral(1,nx),
	 data->Integral(1,nb),
	 datatrue->Integral(1,nb),
	 hxeff->Integral(1,nx)/nx);

  //UnfoldingUtils uu(Adet, data, 0, xini, datatrue, 0);
  UnfoldingUtils uu(Adet, data, 0, xini, 0, hxeff);
  TMatrixD L = uu.LMatrix(nx, UnfoldingUtils::k2DerivBC0);
  uu.SetLMatrix(L);
  uu.SetPrior(xini);
  uu.SetVerbosity(0);

  // GSVD analysis of the system -----------------------------
  // -----------------------------------------------------------------
  GSVDResult* gsvd = uu.GSVDAnalysis(L,110,0,0,"~");
  uu.DrawGSVDPlot(gsvd, 1e-4, 2e2, "~");


  // General-form Tikhonov algorithm using GSVD ----------------------
  // -----------------------------------------------------------------
  int nLambda = 80;
  TVectorD regVector(nLambda);  regVector(0) = 0.;
  for (int k=1; k<nLambda; k++) 
    regVector(k) = 0.005*TMath::Exp(0.15*k);
    //    regVector(k) = 4*k+0.5;
  
  UnfoldingResult rg = uu.UnfoldTikhonovGSVD(gsvd, regVector);
  rg.WRegHist->Scale(double(nx)/(xt2-xt1));
  rg.XRegHist->Scale(double(nx)/(xt2-xt1));
  DrawObject(rg.WRegHist, "surf", "ru_gsvd_wreg2d"); gPad->SetLogy();
  DrawObject(rg.XRegHist, "surf", "ru_gsvd_xreg2d"); gPad->SetLogy();
  DrawObject(rg.GcvCurve, "alp",  "ru_gsvd_gcv"); gPad->SetLogx();
  SetGraphProps(rg.GcvCurve,kMagenta+2,kNone,kMagenta+2,kFullCircle,0.5);
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{min} = %g at k = %d", 
			      rg.lambdaGcv, rg.kGcv));
  TGraph* ggcv = new TGraph(1); 
  ggcv->SetPoint(0,rg.lambdaGcv,rg.GcvCurve->GetY()[rg.kGcv]);
  SetGraphProps(ggcv,kRed,kNone,kRed,kOpenCircle,2);
  ggcv->SetLineWidth(2);
  ggcv->Draw("psame");
  
  DrawObject(rg.LCurve, "alp", "ru_gsvd_lcurve");
  gPad->SetLogx(); gPad->SetLogy();
  SetGraphProps(rg.LCurve,kBlue,kNone,kBlue,kFullCircle,0.5);

  // Only works for square problems
  // // SVD algorithm ---------------------------------------------------
  // // -----------------------------------------------------------------
  // double lambda = 50.0;
  // hSVD = uu.UnfoldSVD(lambda, gsvdHists, "BC0, ~");

  // Chi squared minimization ----------------------------------------
  // -----------------------------------------------------------------
  TVectorD regWts(20);  regWts(0) = 0.;
  for (int k=1; k<20; k++)
    regWts(k) = 1e-9*TMath::Exp(0.5*k); //TMath::Power(2,k)*1e-8;
  uu.SetRegType(UnfoldingUtils::kTotCurv);
  UnfoldingResult cs = uu.UnfoldChiSqMin(regWts);
  cs.XRegHist->Scale(double(nx)/(xt2-xt1));
  DrawObject(cs.XRegHist,"surf", "ru_cs_xreg2d"); gPad->SetLogy();
  DrawObject(cs.LCurve,"alp", "ru_cs_lcurve");
  SetGraphProps(cs.LCurve,kBlue,kNone,kBlue,kFullCircle,0.5);
  hCh2 = cs.XRegHist->ProjectionX(Form("cs%d",10),10,10);

  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nCG = 20;
  UnfoldingResult cg = uu.UnfoldPCGLS(nCG, L,"~");
  cg.WRegHist->Scale(double(nx)/(xt2-xt1));
  cg.XRegHist->Scale(double(nx)/(xt2-xt1));
  DrawObject(cg.XRegHist,"surf", "ru_cg_xreg2d");
  DrawObject(cg.LCurve,"alp", "ru_cg_lcurve");
  SetGraphProps(cg.LCurve,kBlue,kNone,kBlue,kFullCircle,0.5);
  hCG = cg.XRegHist->ProjectionX(Form("cg%d",5),5,5);
  
  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 10;
  //  TH1D* hX0 = data->Clone("hX0"); // use measurement as prior
  //  hX0->Smooth(5);
  UnfoldingResult rl = uu.UnfoldRichardsonLucy(nIterRL);
  rl.XRegHist->Scale(double(nx)/(xt2-xt1));
  hRL = rl.XRegHist->ProjectionX(Form("rl%d",nIterRL),nIterRL,nIterRL);

  // -----------------------------------------------------------------
  // Draw
  // -----------------------------------------------------------------
  SetHistProps(datatrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(data, kBlue, kNone, kBlue, kOpenSquare, 1.0);
  //  SetHistProps(hSVD, kGreen, kNone, kGreen, kOpenCircle, 1.0);
  SetHistProps(rg.hGcv, kGreen+2, kNone, kGreen+2, kOpenCircle, 1.0);
  SetHistProps(hCG, kRed, kNone, kRed, kOpenCircle, 1.5);
  SetHistProps(hCh2, kMagenta+2, kNone, kMagenta+2, kFullCircle, 1.0);
  SetHistProps(hRL, kRed+2, kNone, kRed+2, kOpenCircle, 1.2);

  // Divide by bin width
  datatrue->Scale(1./datatrue->GetBinWidth(1));
  data->Scale(1./data->GetBinWidth(1));
  bini->Scale(1./bini->GetBinWidth(1));
  xini->Scale(1./xini->GetBinWidth(1));
  rg.hGcv->Scale(1./rg.hGcv->GetBinWidth(1));

  // MC and data histos
  DrawObject(hxeff, "", "ru_eff");
  DrawObject(xini, "", "ru_xini_bini");
  bini->SetLineColor(kBlue);
  bini->Draw("same");
  DrawObject(datatrue, "", "ru_data");
  data->Draw("same");

  // Response matrix, as generated (counts)
  DrawObject(Adet,"colz","ru_A");

  // Probability resp. matrix (normalized to efficiency)
  TH2D* hAprob = uu.GetAProbHist();
  hAprob->SetTitle("Detector response (normalized to efficiency);"
		   "measured;true");
  DrawObject(hAprob,"colz","ru_Ahat");

    // Covariance matrices
  if (0) {
    // // Regularization covariance from SVD method
    // TH2D* hXcov = (TH2D*)gsvdHists->FindObject("hXcov1");
    // DrawObject(hXcov, "colz");

    // Statistical error propagation via 100 MC pseudo experiments.
    // Works for any unfolding method, SVD shown here.
    TH2D* hUcov = uu.UnfoldCovMatrix(100, UnfoldingUtils::kGSVDAlgo, 
				     110., "~");
    DrawObject(hUcov, "colz", "ru_mc_xcov");
    // hXcov->Add(hUcov); // Add independent contributions
    for (int j=1; j<=uu.GetN(); j++)
      rg.hGcv->SetBinError(j, TMath::Sqrt(hUcov->GetBinContent(j,j)));
  }

  // Unfolding results
  DrawObject(datatrue, "ep", "ru_results");
  data->Draw("epsame");
  //  hSVD->Draw("epsame");
  rg.hGcv->Draw("epsame");
  hCh2->Draw("epsame");
  hCG->Draw("epsame");
  hRL->Draw("epsame");

  // Animations
  statObjs->Add(datatrue);
  statObjs->Add(data);

  TGraphTime* aRL = Animation(rl.XRegHist, statObjs, "ep", 100, kRed+2, kFullCircle,1.,
			      -10.,0,10,5000);
  TGraphTime* aCG = Animation(cg.XRegHist, statObjs, "ep", 100, kRed, kFullCircle,1.,
			      -10.,0,10,5000);

  if (1) {
    DrawObject(aRL, "", "ru_anim_rl", cList, 700, 500);
    DrawObject(aCG, "", "ru_anim_cg", cList, 700, 500);
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
  for (Int_t i= 0; i<200000; i++) {
    Double_t xt = R.BreitWigner(0.3, 2.5);
    xini->Fill(xt);
    Double_t x = Reconstruct( xt, R );
    if (x != cutdummy) {
      Adet->Fill(x, xt);
      bini->Fill(x);
    }
  }
  
  // Fill the "data" with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<20000; i++) {
    Double_t xt = R.Gaus(0.0, 2.0);
    datatrue->Fill(xt);
    Double_t x = Reconstruct( xt, R );
    if (x != cutdummy) 
      data->Fill(x);
  }
  
  // Fill the data covariance matrix
  for (int i=1; i<=nb; i++) {
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
