
// Copy or symlink the following files to the current directory:
// They are available at https://github.com/andrewadare/utils.git
#include "MatrixUtils.h"
#include "UtilFns.h"

#include "BayesMCFns.h"
#include "TStopwatch.h"
#include "TestProblems.C"

// Example 6.4.1 from "Fully Bayesian Unfolding" arXiv:1201.4612v4

bool doUniformSampling = false;
bool drawReducedVolume = true;
bool printPDFs = true;

// Number of MC likelihood samples
int nMcmcSamples = 1e6;
int nFlatSamples = 5e7;

// Number of simulated counts & event weighting
int nevts = 1e7;
double evtWeight = 1e-3;

// Number of true and reconstructed bins
const int Nt = 14;
const int Nr = 14;

// Number of bins in marginal posterior distributions
int nMcmcBins = 1000;
int nFlatBins = 20;

// Smearing parameters for test problem
double apar = 0.5;
double bpar = 0.1;

// Regularization strength. 
// Try ~1-2 to reduce posterior variance & increase bin-to-bin smoothness.
double alpha = 0.0;

TGraphErrors *DataPoint(TH1 *hD, TH1 *hp, int t, double y=-1);
TGraphErrors *TruePoint(TH1 *hT, TH1 *hp, int t, double y=-1);
TGraphErrors *MD68Point(TH1 *hp, double y = -1);

TGraphAsymmErrors *ReducedSamplingVolume(TH1D **hmp, TGraphAsymmErrors *old);

void BayesUnfoldingExample641()
{

#ifdef __CINT__ // Avoid CINT badness
  Printf("Please compile this script (root BayesUnfoldingExample641.C+)");
  gSystem->Exec("rm AutoDict*");
  gSystem->Exit(0);
#endif

  if (!gROOT->IsBatch())
  {
    Printf("Several canvases coming...adding -b flag.");
    gROOT->SetBatch();
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2f");

  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));

  TRandom3 ran;
  TStopwatch watch; // Watch starts here. A call to Start() would reset it.

  TObjArray *cList = new TObjArray(); // List of drawn canvases --> PDF file

  // Set up the problem
  double bins[Nt+1] = {0};
  for (int j=0; j<=Nt; j++)
    bins[j] = 500*TMath::Exp(0.15*j);

  TestProblem testprob = AtlasDiJetMass(Nt, Nr, bins, bins,
                                        apar, bpar, nevts, evtWeight);
  TH2D *hM   = testprob.Response;
  TH1D *hT   = testprob.xTruth;
  TH1D *hTmc = testprob.xTruthEst;
  TH1D *hMt  = testprob.xIni;
  TH1D *hD   = testprob.bNoisy;
  TH1D *heff = testprob.eff;
  SetHistProps(hT,kRed+2,kNone,kRed+2);
  SetHistProps(hTmc,kRed,kNone,kRed);
  SetHistProps(hD,kBlack,kNone,kBlack,kFullCircle,1.5);

  TMatrixD M   = MatrixUtils::Hist2Matrix(hM);
  TVectorD T   = MatrixUtils::Hist2Vec(hT);   // \hat{T}
  TVectorD Tmc = MatrixUtils::Hist2Vec(hTmc); // \tilde{T}
  TVectorD D   = MatrixUtils::Hist2Vec(hD);
  TVectorD eff = MatrixUtils::Hist2Vec(heff);
  TVectorD Pt  = MatrixUtils::ElemDiv(MatrixUtils::Hist2Vec(hMt), eff); // P(t)
  TMatrixD Prt = MatrixUtils::DivRowsByVector(M, Pt);  // P(r|t)

  // Compute initial sampling volume and do MCMC sampling
  TGraphAsymmErrors *box = HyperBox(hTmc);
  SetGraphProps(box, kGreen+2, kNone, kSpring, kFullSquare, 1.0);

  // Likelihood functor
  LogPoissonLikeFn llfunc(Prt, D);

  // Curvature regularization.
  // Note that this is not directly penalizing curvature of the solution.
  // Instead it smooths the solution divided by the trial spectrum.
  std::vector<double> regpars;
  regpars.push_back(alpha);  // Regularization strength
  for (int i=0; i<box->GetN(); i++)
    regpars.push_back(box->GetY()[i]);

  CurvatureRegFn regfunc(regpars);

  TTree *tmcmc = SampleMH(nMcmcSamples, 1e4, 0.01, box, llfunc, regfunc);

  // Create marginal prob. distributions from MCMC
  std::cout << Form("Marginalizing parameters from Markov chain...")
            << std::flush;

  TH1D *hMCMC[Nt];
  for (int t=0; t<Nt; t++)
  {
    double tlo = box->GetY()[t] - box->GetEYlow()[t];
    double thi = box->GetY()[t] + box->GetEYhigh()[t];
    hMCMC[t] = new TH1D(Form("hMCMC%d",t),"",nMcmcBins, tlo, thi);
    hMCMC[t]->SetTitle(Form("MCMC - point %d;"
                            "entries;"
                            "Marginal posterior probability",t));

    // Marginalize with unit weight when using MCMC, weight by
    // likelihood if sampling was uniform.
    tmcmc->Draw(Form("T%d >> hMCMC%d",t,t), "", "goff");
    hMCMC[t]->Scale(1./hMCMC[t]->Integral(1, nMcmcBins));
    SetHistProps(hMCMC[t], kBlack, kYellow, kBlack, kFullCircle, 1.0);
    hMCMC[t]->GetYaxis()->SetTitleOffset(1.5);
  }
  Printf("Done marginalizing MCMC.");

  // Now compute reduced sampling volume, and do uniform sampling
  TGraphAsymmErrors *rbox = ReducedSamplingVolume(hMCMC, box);
  SetGraphProps(rbox, kBlack, kNone, kNone, kFullSquare, 1.0);
  TH1D *hFlat[Nt];
  if (doUniformSampling)
  {
    TTree *tflat = SampleUniform(nFlatSamples, D, Prt, rbox);
    std::cout << Form("Marginalizing parameters from uniform volume...")
              << std::flush;

    for (int t=0; t<Nt; t++)
    {
      double tlo = rbox->GetY()[t] - rbox->GetEYlow()[t];
      double thi = rbox->GetY()[t] + rbox->GetEYhigh()[t];
      hFlat[t] = new TH1D(Form("hFlat%d",t),"",nFlatBins, tlo, thi);
      hFlat[t]->SetTitle(Form("Uniform sampling - point %d;"
                              "dijet mass (GeV/c^{2});"
                              "Marginal posterior probability",t));

      tflat->Draw(Form("T%d >> hFlat%d",t,t), "L", "goff");
      hFlat[t]->Scale(1./hFlat[t]->Integral(1,nFlatBins));
      SetHistProps(hFlat[t], kBlack, kOrange, kBlack, kFullCircle, 1.0);
    }
    Printf("Done marginalizing uniform volume.");
  }

  // Unfolded spectrum from MCMC
  TGraphErrors *unf1 = new TGraphErrors();
  SetGraphProps(unf1, kBlue, kNone, kBlue, kOpenSquare, 1.5);
  unf1->SetLineWidth(2);
  for (int t=0; t<Nt; t++)
  {
    MaxDensityInterval mdi = GetMDI(hMCMC[t], 0.68);
    unf1->SetPoint(t, hD->GetBinCenter(t+1), mdi.u);
    unf1->SetPointError(t, 0.48*hD->GetBinWidth(t+1), mdi.du);
  }

  // Unfolded spectrum from uniform sampling after volume reduction
  TGraphErrors *unf2 = 0;
  if (doUniformSampling)
  {
    unf2 = new TGraphErrors();
    SetGraphProps(unf2, kRed, kNone, kRed, kOpenSquare, 1.5);
    unf2->SetLineWidth(2);

    for (int t=0; t<Nt; t++)
    {
      MaxDensityInterval mdi = GetMDI(hFlat[t], 0.68);
      unf2->SetPoint(t, hD->GetBinCenter(t+1), mdi.u);
      unf2->SetPointError(t, 0.47*hD->GetBinWidth(t+1), mdi.du);
    }
  }

  Printf("Drawing results...");
  DrawObject(hM, "colz", "matrix", cList, 550, 500);
  gPad->SetLogx();  gPad->SetLogy();  gPad->SetLogz();
  gPad->SetRightMargin(0.15);

  DrawObject(heff, "", "efficiency", cList);

  // Draw marginal dists. from MCMC
  for (int t=0; t<Nt; t++)
  {
    DrawObject(hMCMC[t], "", Form("post_%d", t), cList);
    gPad->SetLeftMargin(0.15);

    if (doUniformSampling)
    {
      hFlat[t]->Scale(1./hFlat[t]->Integral(1, nFlatBins,"width"));
      hFlat[t]->Draw("same");
    }

    double ymin = hMCMC[t]->GetMinimum();
    double ymax = hMCMC[t]->GetMaximum();
    double yDraw = 0.25*(ymax-ymin);
    DataPoint(hD, hMCMC[t], t, 0.75*yDraw)->Draw("ep same");
    TruePoint(hT, hMCMC[t], t, yDraw)->Draw("p same");
    MD68Point(hMCMC[t], yDraw)->Draw("ep same");
  }

  // Result!
  hT->GetYaxis()->SetRangeUser(0.002, 101*nevts*evtWeight);
  DrawObject(hT, "ep", "result", cList);
  gPad->SetLogy();
  box->Draw("e5 same");
  if (doUniformSampling || drawReducedVolume)
    rbox->Draw("e5 same");
  hT->Draw("same");
  hD->Draw("ep same");
  unf1->Draw("ep same");
  if (unf2)
    unf2->Draw("ep same");

  if (printPDFs)
  {
    PrintPDFs(cList, "pdfs"); // Print individuals into ./pdfs dir
    PrintPDF(cList, "pdfs/mcmc_unfold_example"); // Multipage PDF
  }

  Printf("All done.");
  watch.Stop();
  watch.Print();

  return;
}

TGraphAsymmErrors *ReducedSamplingVolume(TH1D **hmp, TGraphAsymmErrors *old)
{
  TGraphAsymmErrors *g = (TGraphAsymmErrors *)old->Clone();

  for (int t=0; t<g->GetN(); t++)
  {

    if (!hmp[t])
      Error("","!hmp[%d]",t);

    MaxDensityInterval mdi = GetMDI(hmp[t], 0.99);
    g->SetPoint(t, g->GetX()[t], mdi.u);
    g->SetPointEYlow(t, mdi.du);
    g->SetPointEYhigh(t, mdi.du);
  }
  return g;
}

TGraphErrors *DataPoint(TH1 *hD, TH1 *hp, int t, double y)
{
  double ymin = hp->GetMinimum();
  double ymax = hp->GetMaximum();
  double yDrawData = y>0? y : ymin + 0.2*(ymax-ymin);
  TGraphErrors *dataPoint = new TGraphErrors();
  dataPoint->SetPoint(0, hD->GetBinContent(t+1), yDrawData);
  dataPoint->SetPointError(0, hD->GetBinError(t+1), 0.0);
  SetGraphProps(dataPoint, kBlack, kNone, kBlack, kFullCircle, 1.0);
  return dataPoint;
}

TGraphErrors *TruePoint(TH1 *hT, TH1 *hp, int t, double y)
{
  double ymin = hp->GetMinimum();
  double ymax = hp->GetMaximum();
  double yDrawData = y>0? y : ymin + 0.2*(ymax-ymin);
  TGraphErrors *truePoint = new TGraphErrors();
  truePoint->SetPoint(0, hT->GetBinContent(t+1), yDrawData);
  SetGraphProps(truePoint, kRed+2, kGray, kRed+2, kOpenCircle, 1.5);
  return truePoint;
}

TGraphErrors *MD68Point(TH1 *hp, double y)
{
  double ymin = hp->GetMinimum();
  double ymax = hp->GetMaximum();
  double yDraw = y>0? y : ymin + 0.25*(ymax-ymin);
  TGraphErrors *g = new TGraphErrors();
  MaxDensityInterval mdi = GetMDI(hp, 0.68);
  g->SetPoint(0, mdi.u, yDraw);
  g->SetPointError(0, mdi.du, 0);
  SetGraphProps(g, kAzure, kAzure, kAzure, kOpenSquare, 1.5);
  return g;
}
