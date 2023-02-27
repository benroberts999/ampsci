#pragma once
#include "qip/Maths.hpp"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <vector>
//
#include <gsl/gsl_errno.h>
#include <iostream>

/*!
@brief Wrappers for returning Spherical Bessel functions.
@details Has an "exact" version, and a faster approx version (good to ~ 1.e-9)
*/
namespace SphericalBessel {

//==============================================================================
template <typename T> T JL(const int L, const T x) {

  // Very low x expansion
  if (std::abs(x) < (T)1.0e-8) {
    if (L == 0)
      return 1.0;
    else
      return 0.0;
  }

  // Low x expansion for first few Ls. Accurate to ~ 10^9
  if (std::fabs(x) < 0.1) {
    if (L == 0)
      return 1.0 - 0.166666667 * (x * x) + 0.00833333333 * qip::pow<4>(x);
    if (L == 1)
      return 0.333333333 * x - 0.0333333333 * (x * x * x) +
             0.00119047619 * qip::pow<5>(x);
    if (L == 2)
      return 0.0666666667 * (x * x) - 0.00476190476 * qip::pow<4>(x) +
             0.000132275132 * qip::pow<6>(x);
    if (L == 3)
      return 0.00952380952 * (x * x * x) - 0.000529100529 * qip::pow<5>(x) +
             0.000012025012 * qip::pow<7>(x);
    if (L == 4)
      return 0.00105820106 * qip::pow<4>(x) - 0.0000481000481 * qip::pow<6>(x) +
             9.25000925e-7 * qip::pow<8>(x);
  }

  // Explicit formalas for first few L's
  if (L == 0)
    return std::sin(x) / x;
  if (L == 1)
    return std::sin(x) / (x * x) - std::cos(x) / x;
  if (L == 2)
    return (3.0 / (x * x) - 1.0) * std::sin(x) / x -
           3.0 * std::cos(x) / (x * x);
  if (L == 3)
    return (15.0 / (x * x * x) - 6.0 / x) * std::sin(x) / x -
           (15.0 / (x * x) - 1.0) * std::cos(x) / x;
  if (L == 4)
    return 5.0 * (-21.0 + 2.0 * x * x) * std::cos(x) / qip::pow<4>(x) +
           (105.0 - 45.0 * x * x + qip::pow<4>(x)) * std::sin(x) /
               (qip::pow<5>(x));

  // If none of above apply, use GSL to calc. accurately
  return (T)gsl_sf_bessel_jl(L, (double)x);
}

//==============================================================================
template <typename T> T exactGSL_JL(int L, T x) {
  // return (T)gsl_sf_bessel_jl(L, (double)x);
  if (L < 0 || x <= 0.0) {
    std::cout << "\n" << L << " " << x << "\n";
  }
  gsl_sf_result result;
  int e = gsl_sf_bessel_jl_e(L, double(x), &result);
  if (e != GSL_SUCCESS) {
    std::cout << "\n" << L << " " << x << "\n";
    return 0.0;
  }
  return result.val;
}

//==============================================================================
//! Creates a vector of Jl(r) for given r
template <typename T>
std::vector<T> fillBesselVec(const int l, const std::vector<T> &xvec) {
  std::vector<T> Jl_vec;
  Jl_vec.reserve(xvec.size());
  for (const auto &x : xvec) {
    Jl_vec.push_back(exactGSL_JL(l, x));
  }
  return Jl_vec;
}

template <typename T>
std::vector<T> fillBesselVec_kr(const int l, const double k,
                                const std::vector<T> &rvec) {
  std::vector<T> Jl_vec;
  Jl_vec.reserve(rvec.size());
  for (const auto &r : rvec) {
    Jl_vec.push_back(exactGSL_JL(l, k * r));
  }
  return Jl_vec;
}

} // namespace SphericalBessel
