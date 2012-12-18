/*
  Code from
  "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  Eur. Phys. J. Plus 127 (2012) 25 
  http://dx.doi.org/10.1140/epjp/i2012-12025-y
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------
  This code is covered by the GNU General Public License:
  http://www.gnu.org/licenses/gpl.html
  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 7 Dec 2011
  $Id$
 */





#include "pValuePoissonError.h"
#include "CompareHistograms.h"

#include<iostream>
using namespace std;

/*
  Given two ROOT histograms (with the same binning!) containing the
  observed and expected counts, create and return a histogram showing
  the significance of their bin-to-bin discrepancies.

  If the histogram representing the expectation (second input
  parameter) has non-zero bin "errors", these are considered the
  standard deviations representing the full uncertainty and the
  significance is computed accordingly, unless this is disabled (third
  parameter).
*/
TH1F* CompareHistograms(TH1* hObs, TH1* hExp,
			bool neglectUncertainty,
			bool variableBinning,
			TH1* hPull)
{
  if (hObs==0 || hExp==0) {
    cerr << "ERROR in CompareHistograms(): invalid input" << endl;
    return 0;
  }
  TString name=hObs->GetName();
  name+="_cmp_";
  name+=hExp->GetName();
  int Nbins = hObs->GetNbinsX();
  if (Nbins != hExp->GetNbinsX()) {
    cerr << "ERROR in CompareHistograms(): different binning" << endl;
    return 0;
  }
  TH1F* hOut = 0;
  if (variableBinning) {
    hOut = new TH1F(name, "",
		    hObs->GetXaxis()->GetNbins(),
		    hObs->GetXaxis()->GetXbins()->GetArray());
  } else {
    hOut = new TH1F(name, "",
		    Nbins,
		    hObs->GetXaxis()->GetXmin(),
		    hObs->GetXaxis()->GetXmax());
  }
  hOut->GetXaxis()->SetTitle( hObs->GetXaxis()->GetTitle() );
  hOut->GetYaxis()->SetTitle("significance");
  hOut->SetFillColor(2);

  for (int i=1; i<=Nbins; ++i) { // SKIP UNDER- AND OVER-FLOWS
    // if (hObs->GetBinContent(i)<=0) continue;

    unsigned nObs = (int) hObs->GetBinContent(i);
    float nExp = hExp->GetBinContent(i);
    float vrnc = hExp->GetBinError(i);
    vrnc *= vrnc; // variance
    float zValue = 0;
    float pValue = 1;
    if (vrnc>0 && !neglectUncertainty) {
      // account for systematic uncertainty
      pValue = pValuePoissonError(nObs, nExp, vrnc);
    } else {
      // assume perfect knowledge of Poisson parameter
      pValue = pValuePoisson(nObs,nExp);
    }
    zValue = pValueToSignificance(pValue, (nObs>nExp));
    if (pValue<0.5) hOut->SetBinContent(i, zValue);
    if (hPull) hPull->Fill(zValue);
  }

  return hOut;

}

