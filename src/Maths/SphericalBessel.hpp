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
template <typename T>
T JL(const int L, const T x) {

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

  // Explicit formulas for first few L's
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
template <typename T>
T exactGSL_JL(int L, T x) {
  return (T)gsl_sf_bessel_jl(L, (double)x);
}

//! Used for frequency-dependent Breit - more generic formula that can be used to evaluate spherical Bessel function of first kind for negative orders
template <typename T>
T exactGSL_JL_alt(int L, T x) {

  // Very low x expansion [keeps terms up to order (wr)^8]
  if (std::abs(x) < (T)1.0e-8) {
    if (L == 0)
      return 1.0;
    else
      return 0.0;
  }

  // series expansion for small wr [keeps every term up to order (wr)^8]
  if (std::fabs(x) < 0.5) {
    if (L == 0) {
      return 1.0 - 0.166666667 * (x * x) + 0.00833333333 * qip::pow<4>(x) -
             0.000198412698 * qip::pow<6>(x) +
             0.00000275573192 * qip::pow<8>(x);
    }
    if (L == 1) {
      return 0.333333333 * x - 0.0333333333 * (x * x * x) +
             0.00119047619 * qip::pow<5>(x) - 0.0000220458553 * qip::pow<7>(x);
    }
    if (L == 2) {
      return 0.0666666667 * (x * x) - 0.00476190476 * qip::pow<4>(x) +
             0.000132275132 * qip::pow<6>(x) -
             0.00000200416867 * qip::pow<8>(x);
    }
    if (L == 3) {
      return 0.00952380952 * (x * x * x) - 0.000529100529 * qip::pow<5>(x) +
             0.0000120250120 * qip::pow<7>(x);
    }
    if (L == 4) {
      return 0.00105820105 * qip::pow<4>(x) - 0.0000481000481 * qip::pow<6>(x) +
             0.000000925000925 * qip::pow<8>(x);
    }
    if (L == 5) {
      return 0.0000962000962 * qip::pow<5>(x) -
             0.0000037000037 * qip::pow<7>(x);
    }
    if (L == 6) {
      return 0.0000074000074 * qip::pow<6>(x) -
             0.00000024666913 * qip::pow<8>(x);
    }
    if (L == 7) {
      return 0.000000493333826 * qip::pow<7>(x);
    }
  }

  // Explicit formulas for first few L's
  if (L == 0) {
    return std::sin(x) / x;
  }
  if (L == 1) {
    return std::sin(x) / (x * x) - std::cos(x) / x;
  }
  if (L == 2) {
    return (3.0 / (x * x) - 1.0) * std::sin(x) / x -
           3.0 * std::cos(x) / (x * x);
  }
  if (L == 3) {
    return (15.0 / (x * x * x) - 6.0 / x) * std::sin(x) / x -
           (15.0 / (x * x) - 1.0) * std::cos(x) / x;
  }
  if (L == 4) {
    return 5.0 * (-21.0 + 2.0 * x * x) * std::cos(x) / qip::pow<4>(x) +
           (105.0 - 45.0 * x * x + qip::pow<4>(x)) * std::sin(x) /
               (qip::pow<5>(x));
  }
  if (L == 5) {
    return (105.0 * x * x - 945.0 - qip::pow<4>(x)) * std::cos(x) /
               qip::pow<5>(x) +
           15.0 * (63.0 - 28.0 * x * x + qip::pow<4>(x)) * std::sin(x) /
               qip::pow<6>(x);
  }
  if (L == 6) {
    return (10395.0 - qip::pow<6>(x) - 4725.0 * x * x +
            210.0 * qip::pow<4>(x)) *
               std::sin(x) / qip::pow<7>(x) +
           21.0 * (60.0 * x * x - 495.0 - qip::pow<4>(x)) * std::cos(x) /
               qip::pow<6>(x);
  }
  if (L == 7) {
    return 7.0 *
               (19305.0 - 4.0 * qip::pow<6>(x) - 8910.0 * x * x +
                450.0 * qip::pow<4>(x)) *
               std::sin(x) / qip::pow<8>(x) +
           (qip::pow<6>(x) + 17325.0 * x * x - 378.0 * qip::pow<4>(x) -
            135135.0) *
               std::cos(x) / qip::pow<7>(x);
  }

  // if (L < 0) {
  //   return 0.0;
  // }

  // return 0.0;

  return (T)(sqrt(M_PI / (2.0 * x)) *
             gsl_sf_bessel_Jnu(L + 1.0 / 2, (double)x));
  //JL(L, x);
}

//! Spherical Bessel functions of second kind
template <typename T>
T exactGSL_YL(int L, T x) {
  return (T)gsl_sf_bessel_yl(L, (double)x);
}

//! Used for frequency-dependent Breit - more generic formula that can be used to evaluate spherical Bessel function of second kind for negative orders
template <typename T>
T exactGSL_YL_alt(int L, T x) {

  // series expansion for small wr [keeps every term up to order (wr)^8]
  if (std::fabs(x) < 0.1) {
    if (L == 0) {
      return -1.0 / x + 0.5 * x - 0.0416666667 * (x * x * x) +
             0.00138888888 * qip::pow<5>(x) - 0.0000248015873 * qip::pow<7>(x);
    }
    if (L == 1) {
      return -1.0 / (x * x) - 0.5 + 0.125 * (x * x) -
             0.00694444444 * qip::pow<4>(x) + 0.000173611111 * qip::pow<6>(x) -
             0.00000248015873 * qip::pow<8>(x);
    }
    if (L == 2) {
      return -3.0 / (x * x * x) - 0.5 / x - 0.125 * x +
             0.0208333333 * (x * x * x) - 0.000868055555 * qip::pow<5>(x) +
             0.0000173611111 * qip::pow<7>(x);
    }
    if (L == 3) {
      return -15.0 / (qip::pow<4>(x)) - 1.5 / (x * x) - 0.125 -
             0.0208333333 * x * x + 0.00260416666 * qip::pow<4>(x) -
             0.0000868055555 * qip::pow<6>(x) +
             0.00000144675925 * qip::pow<8>(x);
    }
    if (L == 4) {
      return -105.0 / (qip::pow<5>(x)) - 7.5 / (x * x * x) - 0.375 / x -
             0.0208333333 * x - 0.00260416666 * (x * x * x) +
             0.000260416666 * qip::pow<5>(x) -
             0.00000723379629 * qip::pow<7>(x);
    }
    if (L == 5) {
      return -945.0 / (qip::pow<6>(x)) - 52.5 / (qip::pow<4>(x)) -
             1.875 / (x * x) - 0.0625 - 0.00260416666 * x * x -
             0.000260416666 * qip::pow<4>(x) +
             0.0000217013888 * qip::pow<6>(x) -
             0.000000516699735 * qip::pow<8>(x);
    }
    if (L == 6) {
      return -10395.0 / (qip::pow<7>(x)) - 472.5 / (qip::pow<5>(x)) -
             13.125 / (x * x * x) - 0.3125 / x - 0.0078125 * x -
             0.000260416666 * (x * x * x) - 0.0000217013888 * qip::pow<5>(x) +
             0.0000015500992 * qip::pow<7>(x);
    }
    if (L == 7) {
      -135135.0 / qip::pow<8>(x) - 5197.5 / qip::pow<6>(x) -
          118.125 / qip::pow<4>(x) - 2.1875 / (x * x) - 0.0390625 -
          0.00078125 * x *x - 0.0000217013888 * qip::pow<4>(x) -
          0.00000155009920 * qip::pow<6>(x) +
          0.0000000968812004 * qip::pow<8>(x);
    }
  }

  // Explicit formulas for first few L's
  if (L == 0) {
    return -std::cos(x) / x;
  }
  if (L == 1) {
    return -std::cos(x) / (x * x) - std::sin(x) / x;
  }
  if (L == 2) {
    return (1.0 - 3.0 / (x * x)) * std::cos(x) / x -
           3.0 * std::sin(x) / (x * x);
  }
  if (L == 3) {
    return (6.0 / x - 15.0 / (x * x * x)) * std::cos(x) / x -
           (15.0 / (x * x) - 1.0) * std::sin(x) / x;
  }
  if (L == 4) {
    return 5.0 * (-21.0 + 2.0 * x * x) * std::sin(x) / qip::pow<4>(x) +
           (45.0 * x * x - 105.0 - qip::pow<4>(x)) * std::cos(x) /
               (qip::pow<5>(x));
  }
  if (L == 5) {
    return 15.0 * (28.0 * x * x - qip::pow<4>(x) - 63.0) * std::cos(x) /
               qip::pow<6>(x) +
           (105.0 * x * x - qip::pow<4>(x) - 945.0) * std::sin(x) /
               qip::pow<5>(x);
  }
  if (L == 6) {
    return (qip::pow<6>(x) - 10395.0 + 4725.0 * x * x -
            210.0 * qip::pow<4>(x)) *
               std::cos(x) / qip::pow<7>(x) +
           21.0 * (60.0 * x * x - 495.0 - qip::pow<4>(x)) * std::sin(x) /
               qip::pow<6>(x);
  }
  if (L == 7) {
    return 7.0 *
               (4.0 * qip::pow<6>(x) + 8910.0 * x * x - 450.0 * qip::pow<4>(x) -
                19305.0) *
               std::cos(x) / qip::pow<8>(x) +
           (qip::pow<6>(x) + 17325.0 * x * x - 378.0 * qip::pow<4>(x) -
            135135.0) *
               std::sin(x) / qip::pow<7>(x);
  }
  // if (L < 0) {
  //   return 0.0;
  // }

  // gsl_set_error_handler_off();
  // return 0.0;

  return (T)(sqrt(M_PI / (2.0 * x)) *
             gsl_sf_bessel_Ynu(L + 1.0 / 2, (double)x));
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

template <typename T>
std::vector<T> fillBesselVec_kr_alt(const int l, const double k,
                                    const std::vector<T> &rvec) {
  std::vector<T> Jl_vec;
  Jl_vec.reserve(rvec.size());
  for (const auto &r : rvec) {
    Jl_vec.push_back(exactGSL_JL_alt(l, k * r));
  }
  return Jl_vec;
}

template <typename T>
std::vector<T> fillSecondBesselVec_kr(const int l, const double w,
                                      const std::vector<T> &rvec) {
  std::vector<T> Yl_vec;
  Yl_vec.reserve(rvec.size());
  for (const auto &r : rvec) {
    Yl_vec.push_back(exactGSL_YL(l, w * r));
  }
  return Yl_vec;
}

template <typename T>
std::vector<T> fillSecondBesselVec_kr_alt(const int l, const double k,
                                          const std::vector<T> &rvec) {
  std::vector<T> Jl_vec;
  Jl_vec.reserve(rvec.size());
  for (const auto &r : rvec) {
    Jl_vec.push_back(exactGSL_YL_alt(l, k * r));
  }
  return Jl_vec;
}

} // namespace SphericalBessel