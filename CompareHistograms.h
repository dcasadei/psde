#ifndef _PSDE_COMPAREHISTOGRAMS_
#define _PSDE_COMPAREHISTOGRAMS_


/*
  Code from 
  "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  Eur. Phys. J. Plus 127 (2012) 25 
  http://dx.doi.org/10.1140/epjp/i2012-12025-y
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 7 Dec 2011
  $Id$
 */


#include "TROOT.h"
#include "TH1.h"



/*
  Given two ROOT histograms (with the same binning!) containing the
  observed and expected counts, create and return a histogram showing
  the significance of their bin-to-bin discrepancies.

  If the histogram representing the expectation (second input
  parameter) has non-zero bin "errors", these are considered the
  standard deviations representing the full uncertainty and the
  significance is computed accordingly, unless this is disabled (third
  parameter).

  The last input pointer is the pull histogram, which is filled with
  all z-values, without discarding the bins for which the p-value is
  larger than 0.5: in case of ideal match, the pulls follow a standard
  normal distribution.  A typical binnin for the pulls is from -5 to
  +5 with 20 bins.
 */
TH1F* CompareHistograms(TH1* hObs=0,  // observed counts
			TH1* hExp=0,  // expectation
			bool neglectUncertainty=false,
			bool variableBinning=false,
			TH1* hPull=0);



#endif
