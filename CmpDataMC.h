#ifndef _PSDE_CMPDATAMC_
#define _PSDE_CMPDATAMC_

/*
  Code from
  "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  Eur. Phys. J. Plus 127 (2012) 25 
  http://dx.doi.org/10.1140/epjp/i2012-12025-y
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 5 Apr 2012
  $ Id: CmpDataMC.h 9 2012-06-24 07:18:53Z casadei $

  -----------------------------------------------------------------

  Use this function when, in addition to the histogram representing
  the observed counts and the histogram representing the expectation,
  a stacked sum of histograms is available which represents all
  contributions to the expectation.

  The significance is computed accounting for the uncertainties
  encoded in the bin "errors" of the expected histogram.
 */


int CmpDataMC(TH1* hObs, // observed counts
	      TH1* hExp, // expectation
	      THStack* stack=0, // (optional) contributions to expectation
	      TLegend* lgstk=0, // (optional) stack legend
	      TString summary="", // multi-page PDF
	      TString eps="",     // single page EPS
	      TString pdf="",     // single page PDF
	      TCanvas* cv=0);

#endif
