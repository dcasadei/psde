/*
 *   Diego Casadei <diego.casadei@cern.ch>
 *   24 Jun 2012
 *   $Id$
 *
 *   ---------------------------------------------------------------
 *
 *   Given the "observed" and "expected" histograms, plot them and
 *   emphasize the significance of the observed deviations from the
 *   expectation.  No "signal" is injected here: the expected
 *   distribution is the true underlying distribution.  The main
 *   purpose is to test the pull distribution.
 *
 *   -----------------------------------------------------------------
 *   This code is covered by the GNU General Public License:
 *   http://www.gnu.org/licenses/gpl.html
 *   -----------------------------------------------------------------
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
 *   [from the command line]$ root -q -b nosignal.C+
 */


#include<iostream>
using namespace std;

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"


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






int nosignal() {

  // histogram containing the expected counts (no uncertainty)
  TH1F* hExp = new TH1F("hExp","Expectation",1000,0,1000);

  // histogram containing the observed counts
  TH1F* hObs = new TH1F("hObs","Observation",1000,0,1000);

  TRandom3 rnd;
  for (int i=0; i<1000; ++i) {
    hExp->SetBinContent(1+i, 10000);
    hObs->SetBinContent(1+i, rnd.Poisson(10000));
  }


  // significance without systematics: 3rd param is "ignore uncertainty"
  TH1F* hPullNoErr = new TH1F("hPullNoErr","Pull distribution;significance",20,-5,5);

  // bin-wise significance
  TH1F* hSigNoErr = CompareHistograms(hObs, hExp, true, false, hPullNoErr);
  hSigNoErr->SetName("hSigNoErr");


  TCanvas* cv = new TCanvas("cv","",600,400);

  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(111);

  cv->SetLogy();
  hPullNoErr->Draw();

  hPullNoErr->Fit("gaus");


  cv->Print("pulls_nosignal.pdf", "pdf");
  cv->Print("pulls_nosignal.eps", "eps");


  return 0;
}


