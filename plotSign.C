/*
 *   Diego Casadei <diego.casadei@cern.ch>
 *   5 Nov 2011
 *   $Id$
 *
 *   After the work by Georgios Choudalakis <georgios.choudalakis@cern.ch>
 *
 *   ---------------------------------------------------------------
 *
 *   Given the "observed" and "expected" histograms, plot them and
 *   emphasize the significance of the observed deviations from the
 *   expectation.
 *
 *   This code is just an example of the use of the functions which
 *   compute the (dis-)agreement between two histograms in units of
 *   significance (number of standard deviations away from the
 *   Gaussian peak).  Look at the functions header files to learn
 *   about their options.
 *
 *   If you use this code, please cite in your bibliography 
 *   "Plotting the Differences Between Data and Expectation"
 *   by Georgios Choudalakis and Diego Casadei
 *   Eur. Phys. J. Plus 127 (2012) 25 
 *   http://dx.doi.org/10.1140/epjp/i2012-12025-y
 *   (http://arxiv.org/abs/1111.2062)
 *
 *   ---------------------------------------------------------------
 *
 *   [from the command line]$ root -q -b plotSign.C+
 */


#include<iostream>
using namespace std;

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"


///
/// Compute the probability of obtaining a deviation as big as the
/// observed one, in the Poisson case of certain or uncertain
/// parameter
///
#include "pValuePoissonError.C"



///
/// Find the significance of the excess/deficit of counts with respect
/// to the expectation.  It returns the histogram of the significance
///
#include "CompareHistograms.C"




int plotSign(TString input="histograms.root") {

  TFile* f = TFile::Open(input);
  if (!f->IsOpen()) {
    cerr << "ERROR: cannot open " << input << endl;
    return -1;
  }

  gStyle->SetOptStat(0);

  // histogram containing the observed counts
  TH1F* hObs = (TH1F*) gDirectory->Get("data");

  // histogram containing the expected counts with their uncertainties
  TH1F* hExp = (TH1F*) gDirectory->Get("bkg");
  hExp->SetMarkerSize(0);

  // significance without systematics: 3rd param is "ignore uncertainty"
  TH1F* hSigNoErr = CompareHistograms(hObs, hExp, true);
  hSigNoErr->SetName("hSigNoErr");

  // significance with systematics: 3rd param is false by default
  TH1F* hSigSyst = CompareHistograms(hObs, hExp);
  hSigSyst->SetName("hSigSyst");

  TCanvas* cv = new TCanvas("cv","",600,600);
  gStyle->SetOptTitle(0);

  // prepare the canvas
  TPad* cv_a = new TPad("cv_a", "",0.0,0.30,1.0,1.0);
  cv_a->SetTopMargin(0.05);
  cv_a->SetBottomMargin(0.001);
  cv_a->Draw();
  TPad* cv_b = new TPad("cv_b", "",0.0,0.0,1.0,0.30);
  cv_b->SetTopMargin(0.0);
  cv_b->SetBottomMargin(0.35);
  cv_b->Draw();

  hExp->SetMarkerSize(0);
  hExp->SetMarkerStyle(0);
  hExp->SetFillColor(kMagenta);

  hSigNoErr->SetLineStyle(1);
  hSigNoErr->SetLineColor(kBlack);
  hSigNoErr->SetFillColor(kRed);
  hSigNoErr->SetMarkerSize(0);
  hSigNoErr->GetYaxis()->SetTitleSize(0.12);
  hSigNoErr->GetYaxis()->SetTitleOffset(0.35);
  hSigNoErr->GetYaxis()->SetLabelSize(0.12);
  hSigNoErr->GetXaxis()->SetLabelSize(0.12);
  hSigNoErr->GetXaxis()->SetTitleSize(0.15);
  hSigNoErr->GetXaxis()->SetTitleOffset(1.0);
  hSigNoErr->GetXaxis()->SetTickLength(0.09);

  TLegend* lg = new TLegend(0.2,0.05,0.7,0.40);
  lg->SetFillStyle(0);
  lg->SetBorderSize(0);
  lg->AddEntry(hObs, "Data", "P");
  lg->AddEntry(hExp, "Expectation", "LF");
  lg->AddEntry((TObject*)0, "", "");
  lg->AddEntry((TObject*)0, "", "");
  lg->AddEntry(hSigNoErr, "Significance (no unc.)", "LF");
  lg->AddEntry(hSigSyst, "Significance with uncertainty", "LF");

  cv_a->cd()->SetLogy();

  TH1F* hExpClone = (TH1F*) hExp->Clone();
  hExpClone->SetFillStyle(0);
  hExpClone->GetYaxis()->SetTitleOffset(0.9);
  hExpClone->Draw("HIST");
  hExp->Draw("E2 SAME");
  hExpClone->Draw("HIST SAME");
  hObs->Draw("SAME");

  lg->Draw();

  cv_b->cd()->SetGridy();
  hSigNoErr->SetAxisRange(-6.6, 6.8, "Y");
  hSigNoErr->Draw("HIST");

  hSigSyst->SetFillColor(kBlue);
  hSigSyst->Draw("HIST,SAME");

  cv->Print("dataVSexpectSyst.pdf", "pdf");
  cv->Print("dataVSexpectSyst.eps", "eps");

  return 0;
}

