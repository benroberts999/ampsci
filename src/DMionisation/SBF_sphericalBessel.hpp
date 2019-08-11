#pragma once
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
/*
Wrappers for returning Spherical Bessel functions.
An "exact" version, and an approx version (good to ~ 1.e-9)
*/

namespace SBF {

//******************************************************************************
template <typename T> T JL(int L, T x) {

  // Very low x expansion
  if (std::fabs(x) < (T)1.e-8) {
    if (L == 0)
      return 1.;
    else
      return 0;
  }

  // Low x expansion for first few Ls. Accurate to ~ 10^9
  if (std::fabs(x) < 0.1) {
    if (L == 0)
      return 1. - 0.166666667 * (x * x) + 0.00833333333 * std::pow(x, 4);
    if (L == 1)
      return 0.333333333 * x - 0.0333333333 * (x * x * x) +
             0.00119047619 * std::pow(x, 5);
    if (L == 2)
      return 0.0666666667 * (x * x) - 0.00476190476 * std::pow(x, 4) +
             0.000132275132 * std::pow(x, 6);
    if (L == 3)
      return 0.00952380952 * (x * x * x) - 0.000529100529 * std::pow(x, 5) +
             0.000012025012 * std::pow(x, 7);
    if (L == 4)
      return 0.00105820106 * std::pow(x, 4) - 0.0000481000481 * std::pow(x, 6) +
             9.25000925e-7 * std::pow(x, 8);
  }

  // Explicit formalas for first few L's
  if (L == 0)
    return sin(x) / x;
  if (L == 1)
    return sin(x) / (x * x) - cos(x) / x;
  if (L == 2)
    return (3. / (x * x) - 1.) * sin(x) / x - 3. * cos(x) / (x * x);
  if (L == 3)
    return (15. / (x * x * x) - 6. / x) * sin(x) / x -
           (15. / (x * x) - 1.) * cos(x) / x;
  if (L == 4)
    return 5. * (-21. + 2. * x * x) * cos(x) / (x * x * x * x) +
           (105. - 45. * x * x + x * x * x * x) * sin(x) / (x * x * x * x * x);

  // If none of above apply, use GSL to calc. accurately
  return gsl_sf_bessel_jl(L, x);
}

//******************************************************************************
template <typename T> T exactGSL_JL(int L, T x) {
  return (T)gsl_sf_bessel_jl(L, (double)x);
}

} // namespace SBF
