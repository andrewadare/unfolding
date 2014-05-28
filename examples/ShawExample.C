#include "UtilFns.h"
#include "MatrixUtils.h"
#include "UnfoldingUtils.h"
#include "TestProblems.h"
#include <TLine.h>

bool printPDFs = true;
const int n=32;                     // number of points

TH1D *hTrue=0;
TH1D *hMeas=0, *hMeasI=0;       // Measured vector b, w/ & w/o noise
TH2D *hResp=0, *hRespT=0;       // Response matrix & transpose
TH1D *hSVD=0;

TObjArray *cList = new TObjArray();
TObjArray *svdResid = new TObjArray();
TObjArray *genSvdAna = new TObjArray();
TObjArray *svdAnim = new TObjArray();
TObjArray *statObjs = new TObjArray();

TLatex lt;
TLatex ltx; // Not NDC

void ShawExample()
{
  
#ifdef __CINT__
  Printf("Please compile this script (root ShawExample.C+) or use ROOT 6.");
  gSystem->Exit(0);
#endif

  gStyle->SetOptStat(0);
  lt.SetNDC();

  // Set the level of Gaussian white noise on b (absolute).  Since
  // cov(err on b) = deltab * I, the prewhitening option "~" is
  // unnecessary for this problem (there is no x^ini).  Only SVD
  // components above this level should contribute to a regularized
  // solution.
  double deltab = 0.05;

  //  TestProblem p = shaw.ShawSystem(n,deltab);

  TestProblem p = ShawSystem(n, deltab);

  hResp  = p.Response;
  hTrue  = p.xTruth;
  hMeas  = p.bNoisy;
  hMeasI = p.bIdeal;

  statObjs->Add(hMeasI);
  statObjs->Add(hMeas);
  statObjs->Add(hTrue);
  SetHistProps(hTrue, kBlack, kNone, kBlack, kFullCircle, 1.0);
  SetHistProps(hMeas, kBlue, kNone, kBlue, kOpenSquare, 1.0);
  SetHistProps(hMeasI, kBlue-1, kNone, kBlue-1, kDot, 1.0);
  hMeasI->SetLineWidth(2);

  UnfoldingUtils uu(hResp, hMeas, 0, 0, hTrue);

  // SVD analysis plot -----------------------------------------------
  // -----------------------------------------------------------------
  // Note that the Fourier coeffs. U'*b level off at deltab.
  SVDResult *svdhists = uu.SVDAnalysis();
  svdhists->U->SetTitle("SVD u-vectors for Shaw test problem;"
                        "Index into x vector;Column index i");
  svdhists->U->GetXaxis()->SetTitleOffset(1.5);
  svdhists->U->GetYaxis()->SetTitleOffset(1.5);
  uu.DrawSVDPlot(svdhists, 1e-5, 1e10);
  DrawObject(svdhists->U, "surf", "shaw_svd_u", cList);

  // GSVD analysis ---------------------------------------------------
  // -----------------------------------------------------------------
  TMatrixD L = uu.LMatrix(n, UnfoldingUtils::k2DerivBC0);
  GSVDResult *gsvd = uu.GSVDAnalysis(L, 0.8);
  gsvd->UHist->SetTitle("GSVD u-vectors for Shaw test problem;"
                        "x;Column index i");
  gsvd->UHist->GetXaxis()->SetTitleOffset(1.5);
  gsvd->UHist->GetYaxis()->SetTitleOffset(1.5);

  DrawObject(gsvd->UHist, "surf", "shaw_gsvd_u", cList);

  uu.DrawGSVDPlot(gsvd, 1e-5, 1e2);
  TLine eLine;
  eLine.DrawLine(0, deltab, n, deltab);
  gPad->SetName("shaw_gsvd_ana");
  cList->Add((TCanvas *)gPad);

  hMeas->SetTitle("Shaw test problem");
  DrawObject(hMeas, "pl", "shaw_test_problem", cList);
  hMeasI->Draw("plsame");
  hTrue->Draw("plsame");
  gsvd->xregHist->Draw("plsame");

  // General-form Tikhonov algorithm using GSVD ----------------------
  // -----------------------------------------------------------------
  int nLambda = 100;
  TVectorD regVector(nLambda);
  for (int k=0; k<nLambda; k++)
    regVector(k) = 0.005*TMath::Exp(0.1*k);

  UnfoldingResult *rg = uu.UnfoldTikhonovGSVD(gsvd, regVector);
  DrawObject(rg->XRegHist, "surf", "shaw_gsvd_solutions", cList);
  gPad->SetLogy();

  SetGraphProps(rg->GcvCurve,kMagenta+2,kNone,kMagenta+2,kFullCircle,0.5);
  DrawObject(rg->GcvCurve, "alp", "shaw_gsvd_gcv", cList);
  gPad->SetLogx();
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{min} = %.2f at k = %d",
                              rg->lambdaGcv, rg->kGcv));
  TGraph *ggcv = new TGraph(1);
  ggcv->SetPoint(0,rg->lambdaGcv,rg->GcvCurve->GetY()[rg->kGcv]);
  SetGraphProps(ggcv,kRed,kNone,kRed,kOpenCircle,2);
  ggcv->SetLineWidth(2);
  ggcv->Draw("psame");

  TGraph *gcvZoom = (TGraph *)rg->GcvCurve->Clone();
  TPad *inset = new TPad("inset", "inset", 0.2, 0.25, 0.6, 0.7, kNone, 1, 0);
  inset->SetFrameBorderMode(0);
  inset->Draw();
  inset->cd();
  TH1F *hfz = inset->DrawFrame(0.005, 0.0092, 1.35, 0.0096);
  hfz->GetXaxis()->SetNdivisions(6,2,0);
  hfz->GetYaxis()->SetNdivisions(5,2,0);
  hfz->GetXaxis()->SetLabelSize(0.08);
  hfz->GetYaxis()->SetLabelSize(0.08);
  inset->SetLeftMargin(0.25);
  inset->SetBottomMargin(0.2);
  inset->Update();

  gcvZoom->Draw("lp same");
  ggcv->Draw("p same");

  SetGraphProps(rg->LCurve,kBlue,kNone,kBlue,kFullCircle,0.5);
  DrawObject(rg->LCurve, "alp", "shaw_gsvd_lcurve", cList, 500, 500);
  rg->LCurve->GetXaxis()->SetRangeUser(0.2, 2);
  rg->LCurve->GetYaxis()->SetRangeUser(0.05, 11);
  gPad->SetLogx();
  gPad->SetLogy();

  TGraph *glc = new TGraph(4);
  glc->SetPoint(0,rg->LCurve->GetX()[rg->kGcv],rg->LCurve->GetY()[rg->kGcv]);
  glc->SetPoint(1,rg->LCurve->GetX()[rg->kRho],rg->LCurve->GetY()[rg->kRho]);
  glc->SetPoint(2,rg->LCurve->GetX()[rg->kStf],rg->LCurve->GetY()[rg->kStf]);
  glc->SetPoint(3,rg->LCurve->GetX()[rg->kLcv],rg->LCurve->GetY()[rg->kLcv]);
  SetGraphProps(glc,kRed,kNone,kRed,kOpenCircle,2);
  glc->SetLineWidth(2);
  glc->Draw("psame");
  ltx.SetTextSize(0.03);
  ltx.DrawLatex(rg->LCurve->GetX()[rg->kGcv],rg->LCurve->GetY()[rg->kGcv],
                Form("1. #lambda_{GCV} = %.2f at k = %d",
                     rg->lambdaGcv, rg->kGcv));
  ltx.DrawLatex(rg->LCurve->GetX()[rg->kLcv],rg->LCurve->GetY()[rg->kLcv],
                Form("2. #lambda_{max} = %.2f at k = %d",
                     rg->lambdaLcv, rg->kLcv));
  ltx.DrawLatex(rg->LCurve->GetX()[rg->kStf],rg->LCurve->GetY()[rg->kStf],
                Form("3. #lambda_{stf} = %.2f at k = %d",
                     rg->lambdaStf, rg->kStf));

  SetGraphProps(rg->LCurvature,kBlue,kNone,kBlue,kFullCircle,0.5);
  DrawObject(rg->LCurvature, "alp", "gsvd_lcurvature", cList);
  gPad->SetLogx();
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{max curvature} = %.2f at k = %d",
                              rg->lambdaLcv, rg->kLcv));

  SetGraphProps(rg->FilterSum,kMagenta+2,kNone,kMagenta+2,kFullCircle,0.5);
  DrawObject(rg->FilterSum, "alp", "shaw_gsvd_fsum", cList);
  gPad->SetLogx();
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{stf} = %.2f at k = %d",
                              rg->lambdaStf, rg->kStf));

  // PCGLS algorithm -------------------------------------------------
  // -----------------------------------------------------------------
  int nIterPCGLS = 7;
  int LMatrixType = UnfoldingUtils::k2DerivBC0;
  UnfoldingResult *cg = uu.UnfoldPCGLS(nIterPCGLS,LMatrixType, "",gsvd);

  SetGraphProps(cg->LCurve, kRed, kNone, kRed, kFullCircle,1.0);
  DrawObject(cg->LCurve, "alp", "shaw_cgls_lcurve", cList);

  ltx.SetTextColor(kRed);
  for (int k=0; k<cg->LCurve->GetN(); k++)
  {
    double x = cg->LCurve->GetX()[k];
    double y = cg->LCurve->GetY()[k];
    ltx.DrawLatex(x, y, Form("%d", k+1));
  }

  SetGraphProps(cg->GcvCurve,kMagenta+1,kNone,kMagenta+1,kFullCircle);
  DrawObject(cg->GcvCurve, "alp", "shaw_cgls_gcv", cList);

  TGraphTime *anim1 = Animation(cg->XRegHist, statObjs, "pl", 200 /*ms*/,
                                kRed, kOpenCircle, 1.0, 0, -1, 1, 4.5);
  DrawObject(anim1, "", "shaw_cgls_anim", 0, 700, 500);
  anim1->SaveAnimatedGif("pdfs/shaw_pcgls_iterations.gif");

  // Richardson-Lucy algorithm ---------------------------------------
  // -----------------------------------------------------------------
  int nIterRL = 100;
  TH1D *prior = 0;
  bool useGsvdSolutionAsPrior = false;
  if (useGsvdSolutionAsPrior)
  {
    prior = (TH1D *)rg->hStf->Clone("rich_lucy_prior");
  }
  else   // Choose a Gaussian prior
  {
    prior = new TH1D("prior", "initial guess",n,0,1);
    for (int i=1; i<=n; i++)
    {
      double x = (i+0.5)*(uu.GetTrueX2()-uu.GetTrueX1())/n;
      double mu = 0.5;
      double sigma = 0.2;
      double norm = 1. / sigma / TMath::Sqrt(TMath::TwoPi());
      prior->SetBinContent(i, norm*TMath::Gaus(x, mu, sigma));
    }
  }

  uu.SetPrior(prior);
  UnfoldingResult *rl = uu.UnfoldRichardsonLucy(nIterRL);
  DrawObject(rl->XRegHist, "surf", "shaw_rl_solutions", cList);

  SetGraphProps(rl->LCurve,kRed+2,kNone,kRed+2,kFullCircle,0.5);
  DrawObject(rl->LCurve, "alp", "shaw_lcurve", cList, 500, 500);
  gPad->SetLogx();
  gPad->SetLogy();

  SetGraphProps(rl->LCurvature,kRed+1,kNone,kRed+1,kFullCircle,0.5);
  DrawObject(rl->LCurvature, "alp", "shaw_rl_lcurvature", cList);
  gPad->SetLogx();
  lt.DrawLatex(0.2, 0.8, Form("#lambda_{max curvature} = %.2f at k = %d",
                              rl->lambdaLcv, rl->kLcv));


  statObjs->Add(prior);
  TGraphTime *anim2 = Animation(rl->XRegHist, statObjs, "pl", 0,
                                kRed+2,kFullCircle, 1, 0, 0, 1, 4.);
  DrawObject(anim2, "", "shaw_rl_anim", 0, 700, 500);


  // Final plotting  -------------------------------------------------
  // -----------------------------------------------------------------
  hResp->SetTitle("Response matrix;x_{meas};x_{true}");
  DrawObject(hResp, "colz", "shaw_response_matrix", cList, 550, 500);
  gPad->SetRightMargin(0.15);

  // Best solutions
  DrawObject(hMeas, "pl", "shaw_solutions_all", cList);
  hMeasI->Draw("plsame");
  hTrue->Draw("plsame");
  rg->hGcv->Draw("psame");
  cg->hLcv = cg->XRegHist->ProjectionX(Form("cg%d",4),4,4);
  cg->hLcv->Draw("psame");
  SetHistProps(rg->hGcv, kGreen+2, kNone, kGreen+2, kFullCircle, 1.0);
  SetHistProps(cg->hLcv, kRed, kNone, kRed, kOpenCircle, 1.0);

  if (printPDFs)
  {
    PrintPDFs(cList, "pdfs");
    PrintPDF(cList, "pdfs/shaw_example");
  }

  return;
}
