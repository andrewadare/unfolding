// BayesMCFns.C
//
// Collection of utility functions implementing the methods described
// in "Fully Bayesian Unfolding" by G. Choudalakis (arXiv:1201.4612v4)
//
// Andrew Adare andrew.adare@colorado.edu
// March 2013

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TTree.h"
#include "TRandom3.h"
#include "ObjectiveFns.h"

#include <iostream>


struct BayesianCredibilityInterval
{
  BayesianCredibilityInterval() :
    u1(0), u2(0), u(0), du(0),
    bin(0), bin1(0), bin2(0),
    probRequested(0), probComputed(0) {}

  BayesianCredibilityInterval(double p) :
    u1(0), u2(0), u(0), du(0),
    bin(0), bin1(0), bin2(0),
    probRequested(p), probComputed(0) {}

  double u1, u2;
  double u;  // Calculated as (u1+u2)/2
  double du; // Calculated as (u2-u1)/2

  int bin, bin1, bin2;
  double probRequested, probComputed;

  TGraphErrors cdf;
};

// Function prototypes
TGraphAsymmErrors *HyperBox(TH1D *h);
TTree *
SampleUniform(int nSamples, TVectorD &D, TMatrixD &Prt,
              TGraphAsymmErrors *box);
TTree *
SampleMH(int nSamples, int nBurnIn, TMatrixD &Prt, TGraphAsymmErrors *box, 
         LogLikeFn &llfunc, LogPrior &priorfunc);
void AssignProposal(const TGraphAsymmErrors *box, const TVectorD &currentvec,
                    TVectorD &newvec);
bool AcceptProposal(double p0, double p1);

void PrintPercentDone(int i, int N, int k);  // Print i/N (in %) every k%.
BayesianCredibilityInterval GetBCI(TH1 *hp, double probFrac);

TGraphAsymmErrors *
HyperBox(TH1D *h)
{
  int Nt = h->GetNbinsX();
  TGraphAsymmErrors *g = new TGraphAsymmErrors(Nt);

  // Hyperbox boundaries
  double min,max,mid;
  for (int t=0; t<Nt; t++)
  {
    mid = h->GetBinContent(t+1);
    min = 1./(t+10) * mid;
    if (min < 0)
      min = 0;
    if (t==0)
      max = 2*mid;
    else
      max = (t+1)*h->GetBinContent(t);

    double ex = h->GetBinWidth(t+1)/2.04;
    g->SetPoint(t,h->GetBinCenter(t+1), mid);
    g->SetPointError(t, ex, ex, mid-min, max-mid);
  }
  return g;
}

TTree *
SampleUniform(int nSamples, TVectorD &D, TMatrixD &Prt, TGraphAsymmErrors *box)
{
  int Nt = box->GetN();
  float Tpoint[Nt], logL=0;
  TRandom3 ran3;
  TVectorD trialT(Nt);
  TVectorD trialR(D.GetNrows());

  TTree *ptree = new TTree("tflat",
                           "posterior probability from uniform sampling");
  for (int t=0; t<Nt; t++)
  {
    ptree->Branch(Form("T%d",t), &Tpoint[t], Form("T%d/F",t));
  }

  ptree->Branch("logL", &logL, "logL/F");

  std::cout << Form("Sampling L(D|T)*pi(T) uniformly...")
            << std::endl;

  for (int i=0; i<nSamples; i++)
  {

    PrintPercentDone(i, nSamples, 1);

    for (int t=0; t<Nt; t++)
    {
      double min = box->GetY()[t] - box->GetEYlow()[t];
      double max = box->GetY()[t] + box->GetEYhigh()[t];
      trialT(t) = ran3.Uniform(min, max);
    }

    trialR = Prt*trialT;
    logL = (float)LogPoissonLikelihood(D,trialR);

    for (int t=0; t<Nt; t++)
    {
      Tpoint[t] = (float)trialT(t);
    }
    ptree->Fill();
  }

  Printf("Filled tree with %lld entries.",ptree->GetEntries());
  return ptree;
}

TTree *
SampleMH(int nSamples, int nBurnIn, TMatrixD &Prt, TGraphAsymmErrors *box, 
         LogLikeFn &llfunc, LogPrior &priorfunc)
{
  // D = measured data, Prt = migration matrix, box = Nt-dim sampling volume.
  // Alpha and Tmc are used to construct the prior for regularization.
  int Nt = box->GetN();
  double p0, p1;      // Current and proposed probabilities
  float Tpoint[Nt], logL;
  TVectorD trialT = MatrixUtils::Graph2Vec(box);
  TVectorD trialR = Prt*trialT;
  TVectorD propT(Nt);

  TTree *ptree = new TTree("tmcmc", "Metropolis-Hastings Markov chain");
  for (int t=0; t<Nt; t++)
  {
    ptree->Branch(Form("T%d",t), &Tpoint[t], Form("T%d/F",t));
  }
  ptree->Branch("logL", &logL, "logL/F");

  // Note that llfunc < 0, priorfunc > 0.
  p0 = llfunc(trialR) - priorfunc(trialT);

  std::cout << Form("Sampling L(D|T)*pi(T) using MCMC...")
            << std::endl;

  for (int i=0; i < nSamples + nBurnIn; i++)
  {
    PrintPercentDone(i, nSamples + nBurnIn, 1);

    // Get a new proposal point (propT). Update trialR.
    AssignProposal(box, trialT, propT);
    trialR = Prt*propT;

    // Note that llfunc < 0, priorfunc > 0.
    p1 = llfunc(trialR) - priorfunc(propT);

    if (AcceptProposal(p0, p1))
    {
      logL = p0 = p1;
      trialT = propT;
      for (int t=0; t<Nt; t++)
      {
        Tpoint[t] = (float)propT(t);
      }
      if (i >= nBurnIn)
        ptree->Fill();
    }
  }

  Printf("Filled tree with %lld entries.", ptree->GetEntries());

  return ptree;
}


void
AssignProposal(const TGraphAsymmErrors *box, const TVectorD &currentvec,
               TVectorD &newvec)
{
  static TRandom3 ran3;

  // Get a proposal point from a small box centered at the current
  // point. Assign the point to newvec.
  for (int t=0; t<box->GetN(); t++)
  {
    double min = box->GetY()[t] - box->GetEYlow()[t];
    double max = box->GetY()[t] + box->GetEYhigh()[t];
    double dT = 0.01*(max - min);

    // Make sure the proposal cell stays within the box, and that it
    // always has the same volume dT^Nt.
    double lo = TMath::Max(min, currentvec(t) - dT/2);
    double hi = TMath::Min(max, currentvec(t) + dT/2);
    if (lo == min)
      hi = lo + dT;
    if (hi == max)
      lo = hi - dT;

    newvec(t) = ran3.Uniform(lo, hi);

    if (newvec(t) < min)
      Warning("","T(%d) = %f < %f", t, newvec(t), min);
    if (newvec(t) > max)
      Warning("","T(%d) = %f > %f", t, newvec(t), max);
  }
  return;
}

bool
AcceptProposal(double p0, double p1)
{
  static TRandom3 ran3;

  if (p1 >= p0)
    return true;

  if (TMath::Log(ran3.Uniform()) < p1-p0)
    return true;

  return false;
}

BayesianCredibilityInterval
GetBCI(TH1 *hp, double probFrac)
{
  // Returns limits of shortest interval in hp containing the
  // probability fraction given by probFrac.
  // The hp histogram is supposed to be a PDF.
  BayesianCredibilityInterval bci(probFrac);

  if (probFrac <= 0. || probFrac >= 1.0)
  {
    Error("BayesianCredibilityInterval()",
          "Requested p %.3f outside 0 < p < 1", probFrac);
    return bci;
  }

  double tot = hp->Integral(1,hp->GetNbinsX());
  if (tot < 0.999 || tot > 1.001)
  {
    Warning("BayesianCredibilityInterval()",
            "PDF histogram integral = %f.\nNormalizing to 1.", tot);
    hp->Scale(1./tot);
  }

  // Create a cumulative probability density graph from hp
  //  TGraph cdf(hp);
  //  int N = cdf.GetN();
  int N = hp->GetNbinsX();
  TGraphErrors cdf(N);
  double psum = 0;
  for (int i=0; i<N; i++)
  {
    //  for (int i=1; i<N; i++) {
    //    cdf.SetPoint(i,cdf.GetX()[i],cdf.GetY()[i]+cdf.GetY()[i-1]);
    psum += hp->GetBinContent(i+1);
    cdf.SetPoint(i, hp->GetBinCenter(i+1), psum);
    cdf.SetPointError(i, 0.5*hp->GetBinWidth(i+1), 0.0);
  }

  // nb bins add up to probFrac, starting at bin i.
  // Initialize to the max. number of bins, then minimize.
  // Assuming bins have uniform width (!)
  int nb = N;

  // Bounds of probFrac starting at bin i
  double p1, p2;

  int i99 = TMath::BinarySearch(N, cdf.GetY(), 0.99);
  bci.u1 = cdf.GetX()[0] - 0.99*cdf.GetEX()[0];
  bci.u2 = cdf.GetX()[i99] + 0.99*cdf.GetEX()[i99];

  // Last bin in probFrac starting at bin i
  int i2 = 0;
  //  Printf("%d", i99);

  for (int i=0; i<i99; i++)
  {
    p1 = cdf.GetY()[i];
    p2 = p1 + probFrac;

    i2 = TMath::BinarySearch(N, cdf.GetY(), p2);

    // for (int j=i+1; j<N; i++) {
    //   if ( cdf.GetY()[j] >= p2 ) {
    //  i2 = j;
    //  break;
    //   }
    // }

    if (i2 > N-2)
      continue;

    //    Printf("%d", i2);

    if (i2-i+1 < nb)
    {
      nb = i2-i+1;
      // bci.u1 = cdf.GetX()[i];
      // bci.u2 = cdf.GetX()[i2];

      bci.u1 = cdf.GetX()[i] - 0.99*cdf.GetEX()[i];
      bci.u2 = cdf.GetX()[i2] + 0.99*cdf.GetEX()[i2];
    }
  }

  bci.u    = (bci.u1+bci.u2)/2;
  bci.du   = (bci.u2-bci.u1)/2;
  bci.bin  = hp->FindBin(bci.u);
  bci.bin1 = hp->FindBin(bci.u1);
  bci.bin2 = hp->FindBin(bci.u2);
  bci.probComputed = hp->Integral(bci.bin1,bci.bin2);
  bci.cdf = cdf;

  double adiff = TMath::Abs(bci.probComputed - bci.probRequested);

  if (adiff > 0.20)
  {
    Error("GetBCI", "Computed (%.2f) vs. requested (%.2f) prob "
          "mismatch for histogram %s.",
          bci.probComputed, bci.probRequested, hp->GetName());
  }
  else if (adiff > 0.05)
  {
    Warning("GetBCI", "Computed probability %.2f for PDF histogram %s "
            "differs from requested %.2f."
            "\nTry narrower bins if they are close.",
            bci.probComputed, hp->GetName(), bci.probRequested);
  }

  return bci;
}

void
PrintPercentDone(int i, int N, int k)
{
  // Print i/N (in %) every k%.

  int percent = i*100./N;
  if (percent % k == 0)
    std::cout << Form("  %d%%\r", percent) << std::flush;
  return;
}
