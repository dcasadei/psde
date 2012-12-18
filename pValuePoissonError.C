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
    return ROOT::Math::inc_gamma(nObs,nExp);
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
    for (unsigned u=1; u<=stop; ++u) {
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
  Normal quantile computed following Peter John Acklam's
  pseudo-code algorithm for rational approximation
  (pjacklam@online.no) 
  http://home.online.no/~pjacklam/notes/invnorm/

  The algorithm below assumes p is the input and x is the output.

   Coefficients in rational approximations.
   a(1) <- -3.969683028665376e+01
   a(2) <-  2.209460984245205e+02
   a(3) <- -2.759285104469687e+02
   a(4) <-  1.383577518672690e+02
   a(5) <- -3.066479806614716e+01
   a(6) <-  2.506628277459239e+00

   b(1) <- -5.447609879822406e+01
   b(2) <-  1.615858368580409e+02
   b(3) <- -1.556989798598866e+02
   b(4) <-  6.680131188771972e+01
   b(5) <- -1.328068155288572e+01

   c(1) <- -7.784894002430293e-03
   c(2) <- -3.223964580411365e-01
   c(3) <- -2.400758277161838e+00
   c(4) <- -2.549732539343734e+00
   c(5) <-  4.374664141464968e+00
   c(6) <-  2.938163982698783e+00

   d(1) <-  7.784695709041462e-03
   d(2) <-  3.224671290700398e-01
   d(3) <-  2.445134137142996e+00
   d(4) <-  3.754408661907416e+00

   Define break-points.

   p_low  <- 0.02425
   p_high <- 1 - p_low

   Rational approximation for lower region.

   if 0 < p < p_low
      q <- sqrt(-2*log(p))
      x <- (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) /
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1)
   endif

   Rational approximation for central region.

   if p_low <= p <= p_high
      q <- p - 0.5
      r <- q*q
      x <- (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q /
           (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1)
   endif

   Rational approximation for upper region.

   if p_high < p < 1
      q <- sqrt(-2*log(1-p))
      x <- -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) /
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1)
   endif

  The absolute value of the relative error (x_approx - x)/x is
  <= 1.15e-9 for all input values p which can be represented in
  IEEE double precision arithmetic.
 */
double pja_normal_quantile(long double p) {

  long double a[6] = {
    -3.969683028665376e+01, // a(1) -> a[0]
     2.209460984245205e+02, // a(2)
    -2.759285104469687e+02, // a(3)
     1.383577518672690e+02, // a(4)
    -3.066479806614716e+01, // a(5)
     2.506628277459239e+00, // a(6) -> a[5]
  };

  long double b[5] = {
    -5.447609879822406e+01, // b(1) -> b[0]
     1.615858368580409e+02, // b(2)
    -1.556989798598866e+02, // b(3)
     6.680131188771972e+01, // b(4)
    -1.328068155288572e+01, // b(5) -> b[4]
  };

  long double c[6] = {
    -7.784894002430293e-03, // c(1) -> c[0]
    -3.223964580411365e-01, // c(2)
    -2.400758277161838e+00, // c(3)
    -2.549732539343734e+00, // c(4)
     4.374664141464968e+00, // c(5)
     2.938163982698783e+00, // c(6) -> c[5]
  };

  long double d[4] = {
     7.784695709041462e-03, // d(1) -> d[0]
     3.224671290700398e-01, // d(2)
     2.445134137142996e+00, // d(3)
     3.754408661907416e+00, // d(4) -> d[3]
  };

  // Define break-points.
  long double p_low  = 0.02425;
  long double p_high = 1 - p_low;

  // output value
  double x=0;

  // Rational approximation for lower region.
  if (0 < p && p < p_low) {
    long double q = sqrt(-2*log(p));
    x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
	((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  // Rational approximation for central region.
  else if (p_low <= p && p <= p_high) {
    long double q = p - 0.5;
    long double r = q*q;
    x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
      (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
  // Rational approximation for upper region.
  else if (p_high < p && p < 1) {
    long double q = sqrt(-2*log(1-p));
    x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }

  return x;
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
  /*
  if (excess) 
    return ROOT::Math::normal_quantile(1-p,1);
  else
    return ROOT::Math::normal_quantile(p,1);
  */

  if (excess) 
    return pja_normal_quantile(1-p);
  else
    return pja_normal_quantile(p);
}

