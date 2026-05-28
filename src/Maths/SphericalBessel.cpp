#include "SphericalBessel.hpp"
#include "qip/Maths.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

namespace SphericalBessel {

//==============================================================================
double jL(int L, double x) {

  if (L == 0 && std::abs(x) < 1.0e-5)
    return 1.0;

  // Low x expansion for first few Ls. Accurate to better than ~ 1e-15
  if (std::abs(x) < 1.0e-5) {
    if (L == 0) {
      return 1.0 - 0.166666666666667 * qip::pow<2>(x);
    } else if (L == 1) {
      return 0.333333333333 * x;
    } else if (L == 2) {
      return 0.06666666667 * qip::pow<2>(x);
    } else if (L == 3) {
      return 0.009523809524 * qip::pow<3>(x);
    } else if (L == 4) {
      return 0.001058201058 * qip::pow<4>(x);
    } else if (L == 5) {
      return 0.00009620009620 * qip::pow<5>(x);
    }
  }

  // Low x expansion for first few Ls. Accurate to better than ~ 1e-15
  if (std::abs(x) < 0.1) {
    if (L == 0) {
      const double x2 = x * x;
      const double x4 = x2 * x2;
      return 1.0 - 0.166666666666667 * x2 + 0.0083333333333333 * x4;
    } else if (L == 1) {
      const double x2 = x * x;
      const double x3 = x2 * x;
      const double x5 = x3 * x2;
      return 0.333333333333 * x - 0.0333333333333 * x3 + 0.001190476190 * x5;
    } else if (L == 2) {
      const double x2 = x * x;
      const double x4 = x2 * x2;
      const double x6 = x4 * x2;
      return 0.06666666667 * x2 - 0.004761904762 * x4 + 0.0001322751323 * x6;
    } else if (L == 3) {
      const double x2 = x * x;
      const double x3 = x2 * x;
      const double x5 = x3 * x2;
      const double x7 = x5 * x2;
      return 0.009523809524 * x3 - 0.0005291005291 * x5 + 0.00001202501203 * x7;
    } else if (L == 4) {
      const double x2 = x * x;
      const double x4 = x2 * x2;
      const double x6 = x4 * x2;
      return 0.001058201058 * x4 - 0.00004810004810 * x6;
    } else if (L == 5) {
      const double x2 = x * x;
      const double x3 = x2 * x;
      const double x5 = x3 * x2;
      const double x7 = x5 * x2;
      return 0.00009620009620 * x5 - 3.700003700e-6 * x7;
    } else if (L == 6) {
      const double x2 = x * x;
      const double x6 = x2 * x2 * x2;
      return 7.400007400e-6 * x6;
    } else if (L == 7) {
      const double x2 = x * x;
      const double x7 = x2 * x2 * x2 * x;
      return 4.933338267e-7 * x7;
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

  if (L < 0)
    return 0.0;

  // If none of above apply, use GSL to calc. accurately
  return gsl_sf_bessel_jl(L, x);
}

//==============================================================================
double exactGSL_JL(int L, double x) { return gsl_sf_bessel_jl(L, x); }

//==============================================================================
double exactGSL_YL(int L, double x) { return gsl_sf_bessel_yl(L, x); }

//==============================================================================
double yL(int L, double x) {

  // series expansion for small wr [keeps every term up to order (wr)^8]
  if (std::abs(x) < 0.1) {
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
      return -135135.0 / qip::pow<8>(x) - 5197.5 / qip::pow<6>(x) -
             118.125 / qip::pow<4>(x) - 2.1875 / (x * x) - 0.0390625 -
             0.00078125 * x * x - 0.0000217013888 * qip::pow<4>(x) -
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

  if (L < 0)
    return 0.0;

  return gsl_sf_bessel_yl(L, x);
}

//==============================================================================

double PhiL(int L, double x, bool tilde) {

  const auto one = tilde ? 0.0 : 1.0;

  if (std::abs(x) < 0.1) {
    if (L == 0) {
      return one - 0.166666667 * x * x + 0.00833333333 * qip::pow<4>(x) -
             0.000198412698 * qip::pow<6>(x) +
             0.00000275573192 * qip::pow<8>(x);
    }
    if (L == 1) {
      return one - 0.1 * x * x + 0.00357142857 * qip::pow<4>(x) -
             0.0000661375661 * qip::pow<6>(x) +
             0.000000751563252 * qip::pow<8>(x);
    }
    if (L == 2) {
      return one - 0.0714285714 * x * x + 0.00198412698 * qip::pow<4>(x) -
             0.0000300625301 * qip::pow<6>(x) +
             0.000000289062789 * qip::pow<8>(x);
    }
    if (L == 3) {
      return one - 0.0555555556 * x * x + 0.00126262626 * qip::pow<4>(x) -
             0.0000161875162 * qip::pow<6>(x) +
             0.000000134895968 * qip::pow<8>(x);
    }
    if (L == 4) {
      return one - 0.0454545455 * x * x + 0.000874125874 * qip::pow<4>(x) -
             0.00000971250971 * qip::pow<6>(x) +
             0.0000000714155126 * qip::pow<8>(x);
    }
    if (L == 5) {
      return one - 0.0384615385 * x * x + 0.000641025641 * qip::pow<4>(x) -
             0.00000628456511 * qip::pow<6>(x) +
             0.0000000413458231 * qip::pow<8>(x);
    }
    if (L == 6) {
      return one - 0.0333333333 * x * x + 0.000490196078 * qip::pow<4>(x) -
             0.0000042999656 * qip::pow<6>(x) +
             0.0000000255950333 * qip::pow<8>(x);
    }
    if (L == 7) {
      return one - 0.0294117647 * x * x + 0.000386996904 * qip::pow<4>(x) -
             0.000003071404 * qip::pow<6>(x) +
             0.000000016692413 * qip::pow<8>(x);
    }
  }

  if (L < 0)
    return 0.0;

  return tilde ?
           (qip::double_factorial(2 * L + 1) / qip::pow(x, L)) * jL(L, x) -
             1.0 :
           (qip::double_factorial(2 * L + 1) / qip::pow(x, L)) * jL(L, x);
}

//==============================================================================
double PsiL(int L, double x, bool tilde) {

  const auto one = tilde ? 0.0 : 1.0;

  if (std::abs(x) < 0.1) {
    if (L == 0) {
      return one - 0.5 * x * x + 0.0416666667 * qip::pow<4>(x) -
             0.00138888889 * qip::pow<6>(x) + 0.0000248015873 * qip::pow<8>(x);
    }
    if (L == 1) {
      return one + 0.5 * x * x - 0.125 * qip::pow<4>(x) +
             0.00694444444 * qip::pow<6>(x) - 0.000173611111 * qip::pow<8>(x);
    }
    if (L == 2) {
      return one + 0.166666667 * x * x + 0.0416666667 * qip::pow<4>(x) -
             0.00694444444 * qip::pow<6>(x) + 0.000289351852 * qip::pow<8>(x);
    }
    if (L == 3) {
      return one + 0.1 * x * x + 0.00833333333 * qip::pow<4>(x) +
             0.00138888889 * qip::pow<6>(x) - 0.000173611111 * qip::pow<8>(x);
    }
    if (L == 4) {
      return one + 0.0714285714 * x * x + 0.00357142857 * qip::pow<4>(x) +
             0.000198412698 * qip::pow<6>(x) + 0.0000248015873 * qip::pow<8>(x);
    }
    if (L == 5) {
      return one + 0.0555555556 * x * x + 0.00198412698 * qip::pow<4>(x) +
             0.0000661375661 * qip::pow<6>(x) +
             0.00000275573192 * qip::pow<8>(x);
    }
    if (L == 6) {
      return one + 0.0454545455 * x * x + 0.00126262626 * qip::pow<4>(x) +
             0.0000300625301 * qip::pow<6>(x) +
             0.000000751563252 * qip::pow<8>(x);
    }
    if (L == 7) {
      return one + 0.0384615385 * x * x + 0.000874125874 * qip::pow<4>(x) +
             0.0000161875162 * qip::pow<6>(x) +
             0.000000289062789 * qip::pow<8>(x);
    }
  }

  if (L < 0)
    return 0.0;

  return tilde ?
           -(qip::pow(x, L + 1) / qip::double_factorial(2 * L - 1)) * yL(L, x) -
             1.0 :
           -(qip::pow(x, L + 1) / qip::double_factorial(2 * L - 1)) * yL(L, x);
}

//==============================================================================
namespace {
// Threshold: if the grid cell at index i has width dx <= lambda/N_per_lambda
// (where lambda = 2 pi / k is the asymptotic wavelength of j_L(kx)),
// the cell resolves the oscillation and we just point-sample j_L. Otherwise
// we cell-average.
constexpr int N_per_lambda = 10;

// Average j_L(k*x) over [x_lo, x_hi] by trapezoidal sub-sampling with n_sub
// intervals. n_sub aims for ~N_per_lambda samples per local wavelength of
// j_L inside the cell, capped at 200. Beyond the cap the cell holds so many
// oscillations that the average is small (Riemann-Lebesgue suppression);
// finer resolution does not change the result meaningfully.
constexpr int N_sub_max = 200;
double cell_average_jL_kx(int l, double k, double x_lo, double x_hi) {
  const double dx = x_hi - x_lo;
  const double samples = double(N_per_lambda) * k * dx / (2.0 * M_PI);
  const int n_sub = std::clamp(int(std::ceil(samples)), 4, N_sub_max);
  const double dxs = dx / double(n_sub);
  double sum = 0.5 * (jL(l, k * x_lo) + jL(l, k * x_hi));
  for (int s = 1; s < n_sub; ++s) {
    sum += jL(l, k * (x_lo + double(s) * dxs));
  }
  return sum / double(n_sub);
}

//==============================================================================
// Core: fill jL_vec[i] with either j_L(k*x[i]) (cell_average == false, or the
// cell at index i is fine enough to resolve j_L) or the cell-averaged value
// (1/dx_i) * integral_{cell_i} j_L(k*x) dx (cell_average == true and the cell
// is too coarse). Cell extents are midpoints to neighbours, asymmetric at the
// endpoints.
void fill_jL_core(int l, double k, const std::vector<double> &xvec,
                  std::vector<double> *jL_vec, bool cell_average) {
  const std::size_t N = xvec.size();
  jL_vec->resize(N);
  if (N == 0)
    return;
  if (N == 1 || !cell_average || k <= 0.0) {
#pragma omp parallel for
    for (std::size_t i = 0; i < N; ++i)
      (*jL_vec)[i] = jL(l, k * xvec[i]);
    return;
  }

  const double dx_target = 2.0 * M_PI / (k * double(N_per_lambda));

#pragma omp parallel for
  for (std::size_t i = 0; i < N; ++i) {
    const double x_lo = (i == 0) ? xvec[0] : 0.5 * (xvec[i - 1] + xvec[i]);
    const double x_hi =
      (i == N - 1) ? xvec.back() : 0.5 * (xvec[i] + xvec[i + 1]);
    const double dx = x_hi - x_lo;
    if (dx <= dx_target) {
      (*jL_vec)[i] = jL(l, k * xvec[i]);
    } else {
      (*jL_vec)[i] = cell_average_jL_kx(l, k, x_lo, x_hi);
    }
  }
}
} // namespace

//==============================================================================
std::vector<double> fillBesselVec(int l, const std::vector<double> &xvec,
                                  bool cell_average) {
  std::vector<double> jl;
  fill_jL_core(l, 1.0, xvec, &jl, cell_average);
  return jl;
}

//==============================================================================
std::vector<double> fillBesselVec_kr(int l, double k,
                                     const std::vector<double> &rvec,
                                     bool cell_average) {
  std::vector<double> jl;
  fill_jL_core(l, k, rvec, &jl, cell_average);
  return jl;
}

//==============================================================================
void fillBesselVec_kr(int l, double k, const std::vector<double> &r,
                      std::vector<double> *jl, bool cell_average) {
  fill_jL_core(l, k, r, jl, cell_average);
}

} // namespace SphericalBessel
