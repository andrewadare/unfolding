TDatime dtime;
TFile* inFileMC    = new TFile("inputs/JetUnfoldPbPbInput.root");
TFile* inFileData  = new TFile("inputs/JetUnfoldPbPbInput.root");
TString outFileName = Form("inputs/andrew-%d.root", dtime.GetDate());
TFile* outFile = new TFile(outFileName, "recreate");

TH1F* hDpt        = (TH1F*)inFileMC->Get("histDeltaPt");
TH2F* hRespDet_in = (TH2F*)inFileMC->Get("histDetectorResponse");
TH1D* hPythia     = hRespDet_in->ProjectionY("hPythia");
TH1D* hMeasIn     = (TH1D*)inFileData->Get("histMeasuredInclusive");

// Set bin width & ranges
const int nMeas = 30;
const int nTrue = 30;
const double ptgen1  = 10;  // generated truth
const double ptgen2  = 160;
const double ptmeas1 = 0;   // measured
const double ptmeas2 = 150;

void CreateResponseMatrix()
{
  Load();
  TH1::AddDirectory(kFALSE); // So histos survive after closing output file
  UnfoldingUtils uu;

  // Create 2D response matrix from dpt distribution. hBkg and
  // hRespDet_in should have identical binning after this step.
  int nm = hRespDet_in->GetNbinsX();
  int nt = hRespDet_in->GetNbinsY();
  TH2D* hBkg = uu.BandedDiagonalMatrix(hDpt, nm, nt);
  hBkg->SetTitle("A_{bkg} from #deltap_{T};"
		 "Reconstructed p_{T} (GeV/c);"
		 "Generated p_{T} (GeV/c)");
  TH2D* hBkgFine = hBkg->Clone("hBkgFine"); // QA
  hBkgFine->SetTitle("A_{bkg} from #deltap_{T} (fine binning)");

  // bkgNorm is projection to true axis. Used later to normalize
  // rebinned hBkg.
  TH1D* bkgNorm = hBkg->ProjectionY("bkgNorm");
  bkgNorm->SetTitle("bkg efficiency");
  TH1D* bkgNormNarrowBins = 
    (TH1D*)bkgNorm->Clone("bkgNormNarrowBins"); // QA

  // Weight hBkg by MC generated spectrum along x/gen/true
  // axis. Accounts for steep spectrum before rebinning.
  for (int j=1; j<=nt; j++) {
    for (int i=1; i<=nm; i++) {
      double mcw = hPythia->GetBinContent(j);
      hBkg->SetBinContent(i,j, mcw*hBkg->GetBinContent(i,j));
    }
  }

  // Rebin bkg TH2 to target bin width
  double rebinMeas = (ptmeas2-ptmeas1)/nMeas/hBkg->GetXaxis()->GetBinWidth(1);
  double rebinTrue = (ptgen2 - ptgen1)/nTrue/hBkg->GetYaxis()->GetBinWidth(1);
  hBkg->Rebin2D(rebinMeas,rebinTrue);
  bkgNorm->Rebin(rebinTrue);
  bkgNorm->Scale(1./rebinTrue);

  // Normalize rebinned hBkg so that each TH2 row j (=TMatrixD cols)
  // sums to content of bkgNorm(j)
  NormalizeXSum(hBkg, bkgNorm);
  TH1D* yProj=(TH1D*)hBkg->ProjectionY("yproj"); // QA

  // Detector response matrix
  TH2D* hDet = (TH2D*)hRespDet_in->Clone("hDet");
  hDet->Rebin2D(rebinMeas,rebinTrue);
  NormalizeXSum(hDet);
  TH1D* DetyProj=(TH1D*)hDet->ProjectionY("detyproj"); // QA
  DetyProj->SetTitle("Detector efficiency");

  // Combine: Afull = Abkg * Adet 
  TH2D* hFull = uu.TH2Product(hBkg, hDet, "hFull");
  hFull->SetTitle("Combined Bkg*Det Response (full range);"
		  "Reconstructed p_{T} (GeV/c);"
		  "Generated p_{T} (GeV/c)");
  
  // Sub-range - Final result!
  int bx1 = hFull->GetXaxis()->FindBin(ptmeas1+0.001);
  int bx2 = hFull->GetXaxis()->FindBin(ptmeas2-0.001);
  int by1 = hFull->GetYaxis()->FindBin(ptgen1+0.001);
  int by2 = hFull->GetYaxis()->FindBin(ptgen2-0.001);
  TH2D* hResp = uu.TH2Sub(hFull, bx1, bx2, by1, by2, "hResp");
  hResp->SetTitle("Bkg*Det Response Matrix;"
		  "Reconstructed p_{T} (GeV/c);"
		  "Generated p_{T} (GeV/c)");
  TH1D* hEffResp=(TH1D*)hResp->ProjectionY("hEffResp"); // QA
  hEffResp->SetTitle("Bkg*Det efficiency");

  // Don't normalize the final hist; eff < 1 should be accounted for
  // during unfolding.  // NormalizeXSum(hResp);

  // Create xini as sub-range of hPythia
  hPythia->SetTitle("MC pt spectrum");

  hPythia->Rebin(rebinTrue);
  TH1D* hXini = new TH1D("hXini", "hXini", nTrue, ptgen1, ptgen2);
  for (int j=1; j<=nTrue; j++) {
    int b = hPythia->FindBin(hXini->GetBinCenter(j));
    double val = hPythia->GetBinContent(b);
    double err = hPythia->GetBinError(b);
    hXini->SetBinContent(j, val);
    hXini->SetBinError(j, err);
  }
  // double areaFrac = hXini->Integral()/hPythia->Integral();
  // hXini->Scale(hPythia->GetEntries()/hXini->GetBinWidth(1)*areaFrac);
  hXini->SetLineColor(kRed);
  hXini->SetLineWidth(2);

  // Rebin measured spectrum
  hMeasIn->Rebin(rebinMeas);
  TH1D* hMeas = new TH1D("hMeas", "Measured pt distribution", nMeas, ptmeas1, ptmeas2);
  for (int i=1; i<=nMeas; i++) {
    int b = hMeasIn->FindBin(hMeas->GetBinCenter(i));
    double val = hMeasIn->GetBinContent(b);
    double err = hMeasIn->GetBinError(b);
    hMeas->SetBinContent(i, val);
    hMeas->SetBinError(i, err);
  }

  DrawObject(hPythia); gPad->SetLogy();
  Printf("hPythia integral %3.2g entries %3.2g", hPythia->Integral(), hPythia->GetEntries());
  hXini->Draw("histsame");
  Printf("hXini integral %3.2g entries %3.2g", hXini->Integral(), hXini->GetEntries());
  DrawObject(hDpt); gPad->SetLogy();
  DrawObject(hBkgFine, "colz");
  DrawObject(hBkg, "colz");
  DrawObject(bkgNorm);
  DrawObject(yProj);
  DrawObject(hDet, "colz");
  Printf("hDet integral %3.2g entries %3.2g", hDet->Integral(), hDet->GetEntries());
  DrawObject(DetyProj);
  DrawObject(hFull, "colz");  
  DrawObject(hResp, "colz");  
  DrawObject(hEffResp);  
  DrawObject(hMeasIn);  
  DrawObject(hMeas);  
  
  hMeasIn->Write();
  hMeas->Write();
  hXini->Write();
  hResp->Write();
  hEffResp->Write();
}

int Load()
{
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->LoadMacro("IOUtilFns.C");
  gROOT->LoadMacro("UnfoldingUtils.C+g");
  return 0;
}

void NormalizeXSum(TH2* hA, TH1* hN=0)
{
  // Normalize x-rows of hA so each row sums to 1.0 (default), or
  // optionally to the value of the jth bin of hN.

  int nx = hA->GetNbinsX();
  int ny = hA->GetNbinsY();

  if (hN)
    if (ny != hN->GetNbinsX())
      Error("NormalizeXSum()", 
	    "ny=%d != %d in hN", nx, hN->GetNbinsX());
  
  // xsum(j) contains sum of x cells in row j
  TVectorD xsum(ny);
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      xsum(j) += hA->GetBinContent(i+1,j+1);
    }
  }

  // Change bin contents of hA to normalized value, which is 1.0 if hN
  // is not passed in.
  for (int j=0; j<ny; j++) {
    double a = hN ? hN->GetBinContent(j+1) : 1.; 
    double w = (xsum(j) > 1e-15) ? a/xsum(j) : a;
    for (int i=0; i<nx; i++) {
      double val = w*hA->GetBinContent(i+1,j+1);
      hA->SetBinContent(i+1,j+1, val);
    }
  }
}
