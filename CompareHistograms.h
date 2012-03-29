#ifndef _COMPAREHISTOGRAMS_
#define _COMPAREHISTOGRAMS_


/*
  Code for "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 7 Dec 2011
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
*/
TH1F* CompareHistograms(TH1* hObs=0, TH1* hExp=0,
			bool neglectUncertainty=false);



#endif
