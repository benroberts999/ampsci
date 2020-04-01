#pragma once
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <vector>

/*!
@brief Wrappers for returning Spherical Bessel functions.
@details Has an "exact" version, and a faster approx version (good to ~ 1.e-9)
*/
namespace SphericalBessel {

//******************************************************************************
template <typename T> T JL(const int L, const T x) {

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
      return 1.0 - 0.166666667 * (x * x) + 0.00833333333 * std::pow(x, 4);
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
    return std::sin(x) / x;
  if (L == 1)
    return std::sin(x) / (x * x) - std::cos(x) / x;
  if (L == 2)
    return (3.0 / (x * x) - 1.) * std::sin(x) / x - 3.0 * std::cos(x) / (x * x);
  if (L == 3)
    return (15.0 / (x * x * x) - 6.0 / x) * std::sin(x) / x -
           (15.0 / (x * x) - 1.) * std::cos(x) / x;
  if (L == 4)
    return 5.0 * (-21.0 + 2.0 * x * x) * std::cos(x) / (x * x * x * x) +
           (105.0 - 45.0 * x * x + x * x * x * x) * std::sin(x) /
               (x * x * x * x * x);

  // If none of above apply, use GSL to calc. accurately
  return (T)gsl_sf_bessel_jl(L, (double)x);
}

//******************************************************************************
template <typename T> T exactGSL_JL(int L, T x) {
  return (T)gsl_sf_bessel_jl(L, (double)x);
}

//******************************************************************************
//! Creates a vector of Jl(r) for given r
template <typename T> //
std::vector<T> fillBesselVec(const int l, const std::vector<T> &xvec) {
  std::vector<T> Jl_vec;
  Jl_vec.reserve(xvec.size());
  for (const auto &x : xvec) {
    Jl_vec.push_back(JL(l, x));
    // Jl_vec.push_back(exactGSL_JL(l, x));
  }
  return Jl_vec;
}

} // namespace SphericalBessel
