#ifndef _PSDE_PVALUEPOISSONERROR_
#define _PSDE_PVALUEPOISSONERROR_

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

  Diego Casadei <casadei@cern.ch> 6 Nov 2011
  $Id$
 */





/*
  p-value for Poisson distribution, no uncertainty on the parameter
*/
double pValuePoisson(unsigned nObs,    // observed counts
		     double nExp);     // Poisson parameter




/*
  p-value for Poisson distribution when there is uncertainty on the
  parameter
*/
double pValuePoissonError(unsigned nObs, // observed counts
			  double E=1,    // expected counts
			  double V=1);   // variance of expectation



/*
  Convert a p-value into a right-tail normal significance, i.e. into
  the number of Gaussian standard deviations which correspond to it.
*/
double pValueToSignificance(double p,          // p-value
			    bool excess=true); // false if deficit


#endif
