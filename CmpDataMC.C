
/*
  Code from
  "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  Eur. Phys. J. Plus 127 (2012) 25 
  http://dx.doi.org/10.1140/epjp/i2012-12025-y
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 5 Apr 2012
  $Id$

  -----------------------------------------------------------------

  Use this function when, in addition to the histogram representing
  the observed counts and the histogram representing the expectation,
  a stacked sum of histograms is available which represents all
  contributions to the expectation.

  The significance is computed accounting for the uncertainties
  encoded in the bin "errors" of the expected histogram.
 */



#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include "CompareHistograms.h"

#include "CmpDataMC.h"

int CmpDataMC(TH1* hObs, // observed counts
	      TH1* hExp, // expectation
	      THStack* stack,  // (optional) contributions to expectation
	      TLegend* legend, // legend
	      TString summary, // multi-page PDF
	      TString eps,     // single page EPS
	      TString pdf,     // single page PDF
	      TCanvas* cv)
{
  if (hObs==0 || hExp==0) {
    cerr << "ERROR in CmpDataMC(): invalid input histograms" << endl;
    return 1;
  }

  if (cv==0) cv = new TCanvas("cv","",600,500);
  cv->Clear();
  cv->Divide();
  gStyle->SetOptStat(0);

  // // significance without systematics: 3rd param is "ignore uncertainty"
  // TH1F* hSigNoErr = CompareHistograms(hObs, hExp, true);
  // hSigNoErr->SetName("hSigNoErr");
  // hSigNoErr->SetLineStyle(1);
  // hSigNoErr->SetLineColor(kBlack);
  // hSigNoErr->SetFillColor(kRed);
  // hSigNoErr->SetMarkerSize(0);
  // hSigNoErr->GetYaxis()->SetTitleSize(0.12);
  // hSigNoErr->GetYaxis()->SetTitleOffset(0.35);
  // hSigNoErr->GetYaxis()->SetLabelSize(0.12);
  // hSigNoErr->GetXaxis()->SetLabelSize(0.12);
  // hSigNoErr->GetXaxis()->SetTitleSize(0.15);
  // hSigNoErr->GetXaxis()->SetTitleOffset(1.0);
  // hSigNoErr->GetXaxis()->SetTickLength(0.09);

  hExp->SetMarkerSize(0);
  hExp->SetMarkerStyle(0);
  hExp->SetFillColor(kCyan-10);

  TH1F* hExpClone = (TH1F*) hExp->Clone();
  hExpClone->SetFillStyle(0);
  hExpClone->GetYaxis()->SetTitleOffset(0.9);

  // significance with systematics: 3rd param is false by default
  TH1F* hSigSyst = CompareHistograms(hObs, hExp);
  if (hSigSyst==0) return 2;
  hSigSyst->SetName("hSigSyst");
  hSigSyst->GetXaxis()->SetLabelSize(0.13);
  hSigSyst->GetYaxis()->SetLabelSize(0.13);
  hSigSyst->GetXaxis()->SetTitleSize(0.15);
  hSigSyst->GetYaxis()->SetTitleSize(0.14);
  hSigSyst->GetXaxis()->SetTitleOffset(1.0);
  hSigSyst->GetYaxis()->SetTitleOffset(0.3);
  hSigSyst->GetXaxis()->SetTickLength(0.09);
  hSigSyst->SetMarkerSize(0);

  // prepare the canvas
  TPad* cv_a = new TPad("cv_a", "",0.0,0.20,1.0,1.0);
  // cv_a->SetTopMargin(0.05);
  // cv_a->SetBottomMargin(0.001);
  cv_a->Draw();
  TPad* cv_b = new TPad("cv_b", "",0.0,0.0,1.0,0.275);
  cv_b->SetTopMargin(0.0);
  cv_b->SetBottomMargin(0.35);
  cv_b->Draw();


  cv_a->cd()->SetLogy();

  double Max=0;
  double Min=1e99;
  if (hExp->GetMaximum()>Max) Max=hExp->GetMaximum();
  if (hObs->GetMaximum()>Max) Max=hObs->GetMaximum();
  if (hExp->GetMinimum()<Min) Min=hExp->GetMinimum();
  if (hObs->GetMinimum()<Min) Min=hObs->GetMinimum();



  if (stack) {
    // stack->SetMinimum(Min);
    stack->SetMaximum(2*Max);
    stack->Draw("HIST");
  } else {
    // hExpClone->SetMinimum(Min);
    hExpClone->SetMaximum(2*Max);
    hExpClone->Draw("HIST");
    hExp->Draw("E2 SAME");
    hExpClone->Draw("HIST SAME");
  }
  hObs->Draw("E SAME");
  if (legend) legend->Draw();


  cv_b->cd()->SetGridy();
  cv_b->cd()->SetGridx();
  // hSigNoErr->SetAxisRange(-6.6, 6.8, "Y");
  // hSigNoErr->Draw("HIST");
  // hSigSyst->SetFillColor(kBlue);
  // hSigSyst->Draw("HIST,SAME");

  // TLegend* lg = new TLegend(0.05,0.05,0.6,0.25);
  // lg->SetFillStyle(0);
  // lg->SetBorderSize(0);
  // lg->AddEntry(hSigNoErr, "No uncertainty", "F");
  // lg->AddEntry(hSigSyst, "With uncertainty", "F");
  // lg->Draw();

  hSigSyst->SetAxisRange(-5.5, 5.5, "Y");
  hSigSyst->Draw("HIST");

  if (! summary.IsNull()) cv->Print(summary, "pdf");
  if (! pdf.IsNull()) cv->Print(pdf, "pdf");
  if (! eps.IsNull()) cv->Print(eps, "eps");

  cv->Clear();
  cv->Divide();
  // delete hExpClone;
  // delete lg;
  // delete cv_a;
  // delete cv_b;

  return 0;
}

