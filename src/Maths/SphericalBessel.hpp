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
//! Spherical Bessel function.
//! For x<0.1 and K<=7, uses expansion accurate to delta~1e-15, eps~1e-3
template <typename T>
T JL(const int L, const T x) {

  // Low x expansion for first few Ls. Accurate to better than ~ 1e-15
  if (std::abs(x) < 0.1) {
    if (L == 0) {
      const T x2 = x * x;
      const T x4 = x2 * x2;
      return 1.0 - 0.166666666666667 * x2 + 0.0083333333333333 * x4;
    } else if (L == 1) {
      const T x2 = x * x;
      const T x3 = x2 * x;
      const T x5 = x3 * x2;
      return 0.333333333333 * x - 0.0333333333333 * x3 + 0.001190476190 * x5;
    } else if (L == 2) {
      const T x2 = x * x;
      const T x4 = x2 * x2;
      const T x6 = x4 * x2;
      return 0.06666666667 * x2 - 0.004761904762 * x4 + 0.0001322751323 * x6;
    } else if (L == 3) {
      const T x2 = x * x;
      const T x3 = x2 * x;
      const T x5 = x3 * x2;
      const T x7 = x5 * x2;
      return 0.009523809524 * x3 - 0.0005291005291 * x5 + 0.00001202501203 * x7;
    } else if (L == 4) {
      const T x2 = x * x;
      const T x4 = x2 * x2;
      const T x6 = x4 * x2;
      return 0.001058201058 * x4 - 0.00004810004810 * x6;
    } else if (L == 5) {
      const T x2 = x * x;
      const T x3 = x2 * x;
      const T x5 = x3 * x2;
      const T x7 = x5 * x2;
      return 0.00009620009620 * x5 - 3.700003700e-6 * x7;
    } else if (L == 6) {
      const T x2 = x * x;
      const T x6 = x2 * x2 * x2;
      return 7.400007400e-6 * x6;
    } else if (L == 7) {
      const T x2 = x * x;
      const T x7 = x2 * x2 * x2 * x;
      return 4.933338267e-7 * x7;
    }
  }

  // If none of above apply, use GSL to calc. accurately
  return (T)gsl_sf_bessel_jl(L, (double)x);
}

//==============================================================================
template <typename T>
T exactGSL_JL(int L, T x) {
  return (T)gsl_sf_bessel_jl(L, (double)x);
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
    Jl_vec.push_back(JL(l, k * r));
  }
  return Jl_vec;
}

template <typename T>
void fillBesselVec_kr(int l, double k, const std::vector<T> &r,
                      std::vector<T> *jl) {
  jl->resize(r.size());
#pragma omp parallel for
  for (std::size_t i = 0; i < r.size(); ++i) {
    (*jl)[i] = JL(l, k * r[i]);
  }
}

} // namespace SphericalBessel
