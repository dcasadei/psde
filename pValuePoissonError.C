/*
  Code for "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 6 Nov 2011
 */




#include<iostream>
#include<cmath>
using namespace std;

#include "TROOT.h"
#include "Math/Math.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"

#include "pValuePoissonError.h"



/*
  p-value for Poisson distribution, no uncertainty on the parameter

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch>   Oct 2011
  Last update: 4 Nov 2011 (using incomplete Gamma from ROOT)

  -----------------------------------------------------------------

  Consider Poi(k|nExp) and compute the p-value which corresponds to
  the observation of nObs counts.

  When nObs > nExp there is an excess of observed events and

    p-value = P(n>=nObs|nExp) = \sum_{n=nObs}^{\infty} Poi(n|nExp)
            = 1 - \sum_{n=0}^{nObs-1} Poi(n|nExp)
            = 1 - e^{-nExp} \sum_{n=0}^{nObs-1} nExp^n / n!

  Otherwise (nObs <= nExp) there is a deficit and

    p-value = P(n<=nObs|nExp) = \sum_{n=0}^{nObs} Poi(n|nExp)
            = e^{-nExp} \sum_{n=0}^{nObs} nExp^n / n!
*/

double pValuePoisson(unsigned nObs,    // observed counts
		     double nExp)      // Poisson parameter
{
  if (nExp==0) return 0.5;
  if (nExp<0) {
    cerr << "ERROR in pValuePoisson(): invalid expectation = " << nExp
	 << " returning 0.5" << endl;
    return 0.5;
  }

  /*
  // Simple recursive formula: Poi(n;nExp) = Poi(n-1;nExp) nExp/n
  // (ROOT independent implementation)
  double p0 = exp(-nExp); // Poi(0;nExp)
  if (nObs>nExp) {// excess
    double pLast = p0;
    double sum = p0;
    for (unsigned k=1; k<=nObs-1; ++k) {
      double p = pLast * nExp / k;
      // cout << Form("Excess: P(%d;%8.5g) = %8.5g and sum = %8.5g",k-1,nExp,pLast,sum) << " -> ";
      sum += p;
      pLast = p;
      // cout << Form("P(%d;%8.5g) = %8.5g and sum = %8.5g",k,nExp,pLast,sum) << endl;
    }
    return 1-sum;
  } else {// deficit
    double pLast = p0;
    double sum = p0;
    for (unsigned k=1; k<=nObs; ++k) {
      // cout << Form("Deficit: P(%d;%8.5g) = %8.5g and sum = %8.5g",k-1,nExp,pLast,sum) << " -> ";
      double p = pLast * nExp / k;
      sum += p;
      pLast = p;
      // cout << Form("P(%d;%8.5g) = %8.5g and sum = %8.5g",k,nExp,pLast,sum) << endl;
    }
    return sum;
  }
  */

  // ROOT provides everything:
  if (nObs>nExp) // excess
    return 1-ROOT::Math::inc_gamma_c(nObs,nExp);
  else // deficit
    return ROOT::Math::inc_gamma_c(nObs+1,nExp);
}





/*

  p-value for Poisson distribution when there is uncertainty on the
  parameter

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch>  6 Nov 2011
  Last update: 3 Mar 2012 (logarithms used only for big numbers)

  -----------------------------------------------------------------

  Consider Poi(k|nExp) and compute the p-value which corresponds to
  the observation of nObs counts, in the case of uncertain nExp whose
  variance is provided.

  The prior for nExp is a Gamma density which matches the expectation
  and variance provided as input.  The marginal model is provided by
  the Poisson-Gamma mixture, which is used to compute the p-value.

  Gamma density: the parameters are
   * a = shape param  [dimensionless]
   * b = rate param   [dimension: inverse of x]

    nExp ~ Ga(x|a,b) = [b^a/Gamma(a)] x^{a-1} exp(-bx)

  One has E[x] = a/b and V[x] = a/b^2 hence
   * b = E/V
   * a = E*b

  The integral of Poi(n|x) Ga(x|a,b) over x gives the (marginal)
  probability of observing n counts as

                b^a [Gamma(n+a) / Gamma(a)]
    P(n|a,b) = -----------------------------
                       n! (1+b)^{n+a}

  When nObs > nExp there is an excess of observed events and

    p-value = P(n>=nObs) = \sum_{n=nObs}^{\infty} P(n)
            = 1 - \sum_{n=0}^{nObs-1} P(n)

  Otherwise (nObs <= nExp) there is a deficit and

    p-value = P(n<=nObs) = \sum_{n=0}^{nObs} P(n)

  To compute the sum, we use the following recurrent relation:

    P(n=0) = [b/(1+b)]^a
    P(n=1) = [b/(1+b)]^a a/(1+b) = P(n=0) a/(1+b)
    P(n=2) = [b/(1+b)]^a a/(1+b) (a+1)/[2(1+b)] = P(n=1) (a+1)/[2(1+b)]
    ...        ...
    P(n=k) = P(n=k-1) (a+k-1) / [k(1+b)]

  and to avoid rounding errors, we work with logarithms.
*/

double pValuePoissonError(unsigned nObs, // observed counts
			  double E,      // expected counts
			  double V)      // variance of expectation
{
  if (E<=0 || V<=0) {
    cerr << "ERROR in pValuePoissonError(): expectation and variance must be positive. "
	 << "Returning 0.5" << endl;
    return 0.5;
  }
  double B = E/V;
  double A = E*B;

  // relative syst = sqrt(V)/E = 1/sqrt(A)
  // relative stat = 1/sqrt(nObs)
  // if syst < 0.1*stat there is no need for syst:
  // save a bit of CPU time :-)
  // if (A>100*nObs) return pValuePoisson(nObs,E); // UNCOMMENT TO SPEED-UP

  if (A>100) { // need to use logarithms

    unsigned stop=nObs;
    if (nObs>E) --stop;

    /// NB: must work in log-scale otherwise troubles!
    double logProb = A*log(B/(1+B));
    double sum=exp(logProb); // P(n=0)
    for (unsigned u=1; u<stop; ++u) {
      logProb += log((A+u-1)/(u*(1+B)));
      sum += exp(logProb);
    }
    if (nObs>E)  // excess
      return 1-sum;
    else  // deficit
      return sum;

  } else {

    // Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
    double p0 = pow(B/(1+B),A); // P(0;A,B)
    double nExp = A/B;
    if (nObs>nExp) {// excess
      double pLast = p0;
      double sum = p0;
      for (unsigned k=1; k<=nObs-1; ++k) {
	double p = pLast * (A+k-1) / (k*(1+B));
	// cout << Form("Excess: P(%d;%8.5g) = %8.5g and sum = %8.5g",k-1,nExp,pLast,sum) << " -> ";
	sum += p;
	pLast = p;
	// cout << Form("P(%d;%8.5g) = %8.5g and sum = %8.5g",k,nExp,pLast,sum) << endl;
      }
      return 1-sum;
    } else {// deficit
      double pLast = p0;
      double sum = p0;
      for (unsigned k=1; k<=nObs; ++k) {
	// cout << Form("Deficit: P(%d;%8.5g) = %8.5g and sum = %8.5g",k-1,nExp,pLast,sum) << " -> ";
	double p = pLast * (A+k-1) / (k*(1+B));
	sum += p;
	pLast = p;
	// cout << Form("P(%d;%8.5g) = %8.5g and sum = %8.5g",k,nExp,pLast,sum) << endl;
      }
      return sum;
    }
  }

}



/*
  Convert a p-value into a right-tail normal significance, i.e. into
  the number of Gaussian standard deviations which correspond to it.

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch>  Oct 2011
*/

#include<iostream>
using namespace std;

#include "TROOT.h"
#include "Math/Math.h"
#include "Math/QuantFuncMathCore.h"

double pValueToSignificance(double p,     // p-value
			    bool excess)  // false if deficit
{
  if (p<0 || p>1) {
    cerr << "ERROR: p-value must belong to [0,1] but input value is " << p << endl;
    return 0;
  }

  if (excess) 
    return ROOT::Math::normal_quantile(1-p,1);
  else
    return ROOT::Math::normal_quantile(p,1);
}

