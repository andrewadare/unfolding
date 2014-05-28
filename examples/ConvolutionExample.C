#include "UtilFns.h"
#include "MatrixUtils.h"
#include "UnfoldingUtils.h"
#include "TestProblems.h"
#include "TLegend.h"

bool printPDFs = true;   // Print canvases to pdf
int m = 100, n = 25;     // Number of measured and true bins
double x1 = 0., x2 = 5.; // Range of solution

// Response matrix, truth model, truth model estimator, and data.
TH2D *hResp=0;
TH1D *hTrue=0;
TH1D *hTrueData=0;
TH1D *hMeas=0;

// Unfolding solutions
// TH1D* hSVD=0;
TH1D *hRL=0;
TH1D *hCh2=0;


TObjArray *cList = new TObjArray(); // Canvas container

void ConvolutionExample()
{
  gStyle->SetBarWidth(0.8);
  gStyle->SetOptStat(0);
  TLatex lt;
  lt.SetNDC();

  // Create test problem from truth and convolver TF1s
  UnfoldingUtils utils;

  TF1 *fTrue = new TF1("fTrue", "TMath::Exp(-x)", x1,x2);
  TF1 *fGaus = new TF1("fGaus", "[0]*TMath::Gaus(x,0,1,1)", -x2, x2);
  fGaus->SetParameter(0,(x2-x1)/m);

  TestProblem t = MonteCarloConvolution(m, n, x1, x2, x1, x2, fTrue, fGaus, 10000);

  hResp = t.Response;
  hTrue = t.xTruth;
  hTrueData = t.xTruthEst;
  hMeas = t.bNoisy;

  // Create UnfoldingUtils instance
  // Use measured data as guess for solution. Clone possible only for m = n.
  TH2D *hMeasCov = 0;
  TH1D *hEff = hResp->ProjectionY("hEff");
  UnfoldingUtils uu(hResp, hMeas, hMeasCov, 0, 0, hEff);

  Printf("||e||_2 = %g, mean = %g, RMS = %g",
         uu.GetbErrNorm(), uu.GetbErrMean(), uu.GetbErrRMS());

  // SVD analysis ----------------------------------------------------
  // -----------------------------------------------------------------
  SVDResult *svd = uu.SVDAnalysis();

  // Draw s.v. coefficients |u_i'*b|
  cList->Add(uu.DrawSVDPlot(svd, 0.2, 1e9));
  gPad->SetName("conv_svd_ana");

  // GSVD analysis ---------------------------------------------------
  // -----------------------------------------------------------------
  TMatrixD L = uu.LMatrix(n, UnfoldingUtils::k2DerivNoBC);
  GSVDResult *gsvd = uu.GSVDAnalysis(L,0.5,0,0,"");
  TString opt = "SPEC a(30,30,180) dm(0,6)";
  gsvd->UHist->SetTitle("GSVD u-vectors for Gaussian convolution matrix;"
                        "x;Column index i");
  gsvd->UHist->GetXaxis()->SetTitleOffset(1.5);
  gsvd->UHist->GetYaxis()->SetTitleOffset(1.5);
  DrawObject(gsvd->UHist, "surf", "conv_gsvd_uvectors", cList);
  cList->Add(uu.DrawGSVDPlot(gsvd, 1e-5, 1e6));
  gPad->SetName("conv_gsvd_ana");

  // General-form Tikhonov algorithm using GSVD ----------------------
  // -----------------------------------------------------------------
  int nLambda = 100;
  TVectorD regVector(nLambda);
  for (int k=0; k<nLambda; k++)
    regVector(k) = 0.01*(k+1);

  UnfoldingResult *rg = uu.UnfoldTikhonovGSVD(gsvd, regVector);
  Info("", "Finished GSVD analysis.");

  DrawObject(rg->XRegHist, "surf"/*"SPEC dm(0,2)"*/, "conv_gsvd_x", cList);

  DrawObject(rg->GcvCurve, "alp", "conv_gsvd_gcv", cList);
  SetGraphProps(rg->GcvCurve,kMagenta+2,kNone,kMagenta+2,kFullCircle,0.5);
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{min} = %g at k = %d",
                              rg->lambdaGcv, rg->kGcv));
  TGraph *ggcv = new TGraph(1);
  ggcv->SetPoint(0,rg->lambdaGcv,rg->GcvCurve->GetY()[rg->kGcv]);
  SetGraphProps(ggcv,kRed,kNone,kRed,kOpenCircle,2);
  ggcv->SetLineWidth(2);
  ggcv->Draw("psame");

  DrawObject(rg->LCurve, "alp", "conv_gsvd_lcurve", cList);
  SetGraphProps(rg->LCurve,kBlue,kNone,kBlue,kFullCircle,0.5);
  TGraph *ggl = new TGraph(1);
  ggl->SetPoint(0,rg->LCurve->GetX()[rg->kGcv],rg->LCurve->GetY()[rg->kGcv]);
  SetGraphProps(ggl,kRed,kNone,kRed,kOpenCircle,2);
  ggl->SetLineWidth(2);
  ggl->Draw("psame");

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 200;
  UnfoldingResult *rl = uu.UnfoldRichardsonLucy(nIterRL);
  DrawObject(rl->XRegHist,"surf");
  DrawObject(rl->LCurve, "alp", "conv_richlucy_lcurve", cList);
  hRL = rl->XRegHist->ProjectionX(Form("rl%d",nIterRL),nIterRL,nIterRL);
  hRL->Scale(1./hRL->GetBinWidth(1));

  // Chi squared minimization ----------------------------------------
  // -----------------------------------------------------------------
  int nl = 50;
  TVectorD regWts(nl);
  for (int k=0; k<nl; k++)
    regWts(k) = (k+1)*1e-5;
  uu.SetRegType(UnfoldingUtils::kTotCurv);
  UnfoldingResult *cs = uu.UnfoldChiSqMin(regWts);
  DrawObject(cs->XRegHist,"surf");
  hCh2 = cs->XRegHist->ProjectionX(Form("cs%d",30),30,30);
  hCh2->Scale(1./hCh2->GetBinWidth(1));
  DrawObject(cs->LCurve,"alp");
  SetGraphProps(cs->LCurve,kBlue,kNone,kBlue,kFullCircle,0.5);

  // -----------------------------------------------------------------
  // Draw
  // -----------------------------------------------------------------

  // Rescale for drawing (do not rescale before unfolding!)
  hTrue->Scale(1./hTrue->GetBinWidth(1));
  hTrueData->Scale(1./hTrue->GetBinWidth(1));
  hMeas->Scale(1./hMeas->GetBinWidth(1));

  SetHistProps(hMeas, kBlue, kNone, kBlue, kOpenSquare, 0.8);
  SetHistProps(hTrue, kBlack, kNone, kBlack, kFullCircle, 0.8);
  SetHistProps(hTrueData, kBlack, kNone, kBlack, kFullCircle, 0.8);
  hTrue->SetLineWidth(2);
  //  SetHistProps(hSVD, kGreen+2, kNone, kGreen+2, kOpenCircle, 1.0);
  SetHistProps(rg->hGcv, kGreen+2, kNone, kGreen+2, kFullSquare, 1.2);
  // SetHistProps(hCG, kRed, kNone, kRed, kOpenCircle, 1.5);
  SetHistProps(hCh2, kMagenta+2, kNone, kMagenta+2, kFullCircle, 1.0);
  SetHistProps(hRL, kRed+2, kNone, kRed+2, kOpenCircle, 1.2);

  // Draw response matrix
  DrawObject(hResp, "col", "conv_response", cList, 500, 500);

  // Reconstruction efficiency
  hEff->SetTitle("Efficiency;t;#epsilon_{j}");
  DrawObject(hEff, "", "conv_resp_eff", cList, 700, 500);

  // Draw the problem without solutions
  DrawObject(hTrue, "l", "conv_problem_setup", cList, 500, 500);
  hTrueData->Draw("epsame");
  hMeas->Draw("epsame");
  TLegend *l0 = new TLegend(0.5, 0.6, 0.99, 0.99,
                            "#splitline{Exponential truth,}{Gaussian convolution}");
  l0->SetFillColor(kNone);
  l0->AddEntry(hTrue, "Data input model", "l");
  l0->AddEntry(hMeas, "Measurement", "epl");
  l0->Draw();

  // Draw true & meas
  DrawObject(hTrue, "l", "conv_problem", cList, 700, 500);
  hTrue->GetYaxis()->SetRangeUser(0., 1.2*hTrue->GetMaximum());
  hMeas->Draw("epsame");
  rg->hGcv->Scale(1./rg->hGcv->GetBinWidth(1));
  rg->hGcv->Draw("psame");
  hRL->Draw("epsame");
  hCh2->Draw("epsame");
  TLegend *l1 = new TLegend(0.5, 0.6, 0.99, 0.99);
  l1->SetFillColor(kNone);
  l1->AddEntry(hTrue, "Data input model", "l");
  l1->AddEntry(hMeas, "Measurement", "epl");
  l1->AddEntry(rg->hGcv, "Tikhonov GSVD", "epl");
  l1->AddEntry(hRL, Form("Richardson-Lucy alg."), "epl");
  l1->AddEntry(hCh2, Form("#chi^{2} minimization method"), "epl");
  l1->Draw();

  // Draw a few of the left singular vectors
  TH1D *hu[999];
  for (int i=0; i<n; i++)
  {
    hu[i] = svd->U->ProjectionX(Form("u_%d",i),i+1,i+1);
    hu[i]->SetLineWidth(2);
    hu[i]->SetTitle(Form("Left singular vector u_{%d};Column index",i));
  }
  SetHistProps(hu[n-1], kGray, kNone, kGray, kFullCircle, 0.8);
  DrawObject(hu[n-1], "pl", "conv_u", cList);

  hu[n-1]->SetTitle(Form("Left Singular Vectors "
                         "#color[4]{u_{1}}, "
                         "#color[2]{u_{2}}, "
                         "#color[418]{u_{3}}, "
                         "#color[922]{u_{%d}}", n));

  SetHistProps(hu[2], kGreen+1, kNone, kGreen+1, kFullCircle, 0.8);
  hu[2]->Draw("plsame");
  SetHistProps(hu[1], kRed, kNone, kRed, kFullCircle, 0.8);
  hu[1]->Draw("plsame");
  SetHistProps(hu[0], kAzure-2, kNone, kAzure-2, kFullCircle, 0.8);
  hu[0]->Draw("plsame");

  // Draw the s.v. spectrum
  DrawObject(svd->sigma, "p", "conv_sigma", cList, 500, 500);
  gPad->SetLogy();

  if (printPDFs)
  {
    PrintPDFs(cList, "pdfs");
    PrintPDF(cList, "pdfs/conv_example");
  }

  return;
}
