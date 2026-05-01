#pragma once
#include "LinAlg/Matrix.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
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

  return (T)(sqrt(M_PI / (2.0 * x)) *
             gsl_sf_bessel_Jnu(L + 1.0 / 2, (double)x));
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
        0.00000155009920 * qip::pow<6>(x) + 0.0000000968812004 * qip::pow<8>(x);
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
    // Jl_vec.push_back(exactGSL_JL(l, x));
    Jl_vec.push_back(JL(l, x));
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

//==============================================================================
//! Spherical Bessel Lookup table; in the form j_L(qr) = J[L][q][r].
/*!
@details
 - Allows 'first match', 'nearest', and 'interpolation' lookup
 - All of these are only accurate if grid dense enough
 - Ideally, use exact same grid to extract values as used to create them
*/
class JL_table {
public:
  JL_table() {}
  JL_table(int max_L, const std::vector<double> &q,
           const std::vector<double> &r) {
    fill(max_L, q, r);
  }

  //! If derivatives required for degree L, max_L = L+1
  void fill(int max_L, const std::vector<double> &q,
            const std::vector<double> &r) {
    m_q = q;
    m_J_L_q.resize(std::size_t(max_L + 1), q.size());
    m_J_L_q_on_qr.resize(std::size_t(max_L + 1), q.size());
    using namespace qip::overloads;
#pragma omp parallel for collapse(2)
    for (auto L = 0ul; L <= std::size_t(max_L); ++L) {
      for (auto iq = 0ul; iq < m_J_L_q.cols(); ++iq) {
        const auto tq = q[iq];
        // Fill jL(qr)
        m_J_L_q[L][iq] = fillBesselVec_kr(int(L), tq, r);
        // Store jL(qr)/qr
        m_J_L_q_on_qr[L][iq] = m_J_L_q[L][iq] / (tq * r);
      }
    }
  }

  //----------------------------------------------------------------------------
  //! Direct access to jL(q*r) for specific q grid index
  const std::vector<double> &at(std::size_t L, std::size_t iq) const {
    return m_J_L_q.atc(L, iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q_i*r) for the first grid point q_i such that q_i >= q
  const std::vector<double> &jL(int L, double q) const {
    // The '-1' here ensures we always get a valid index, even if q > m_q.back()
    const auto it = std::lower_bound(m_q.begin(), m_q.end() - 1, q);
    const auto iq = std::size_t(std::distance(m_q.begin(), it));
    return m_J_L_q.atc(std::size_t(L), iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q_i*r) for the grid point q_i nearest to the requested q.
  const std::vector<double> &jL_nearest(int L, double q) const {

    // Clamp
    if (q <= m_q.front())
      return m_J_L_q.atc(std::size_t(L), 0);
    if (q >= m_q.back())
      return m_J_L_q.atc(std::size_t(L), m_q.size() - 1);

    // Find i such that m_q[i] <= q < m_q[i+1]
    const auto it = std::lower_bound(m_q.begin() + 1, m_q.end() - 1, q);
    const auto i = std::size_t(std::distance(m_q.begin(), it)) - 1;

    // Pick nearest
    const auto iq = (q - m_q[i] < m_q[i + 1] - q) ? i : i + 1;

    return m_J_L_q.atc(std::size_t(L), iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q_i*r)/(q_i*r) for the grid point q_i nearest to the requested q.
  const std::vector<double> &jL_on_qr_nearest(int L, double q) const {

    // Clamp
    if (q <= m_q.front())
      return m_J_L_q_on_qr.atc(std::size_t(L), 0);
    if (q >= m_q.back())
      return m_J_L_q_on_qr.atc(std::size_t(L), m_q.size() - 1);

    // Find i such that m_q[i] <= q < m_q[i+1]
    const auto it = std::lower_bound(m_q.begin() + 1, m_q.end() - 1, q);
    const auto i = std::size_t(std::distance(m_q.begin(), it)) - 1;

    // Pick nearest
    const auto iq = (q - m_q[i] < m_q[i + 1] - q) ? i : i + 1;

    return m_J_L_q_on_qr.atc(std::size_t(L), iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q*r) interpolated linearly between grid points.
  //! nb: assumes q grid dense enough that linear interp is reasonable!
  std::vector<double> jL_interp(int L, double q) const {

    // Clamp q to the tabulated range
    if (q <= m_q.front())
      return m_J_L_q.atc(std::size_t(L), 0);
    if (q >= m_q.back())
      return m_J_L_q.atc(std::size_t(L), m_q.size() - 1);

    // Find i such that m_q[i] <= q < m_q[i+1]
    const auto it = std::lower_bound(m_q.begin() + 1, m_q.end() - 1, q);
    const auto i = std::size_t(std::distance(m_q.begin(), it)) - 1;

    // Linear interpolation weight
    // const double t = (q - m_q[i]) / (m_q[i + 1] - m_q[i]);

    const double logq = std::log10(q);
    const double logq0 = std::log10(m_q[i]);
    const double logq1 = std::log10(m_q[i + 1]);
    const double t = (logq - logq0) / (logq1 - logq0);

    const auto &v0 = m_J_L_q.atc(std::size_t(L), i);
    const auto &v1 = m_J_L_q.atc(std::size_t(L), i + 1);

    using namespace qip::overloads;
    return (1.0 - t) * v0 + t * v1;
  }

  //----------------------------------------------------------------------------
private:
  LinAlg::Matrix<std::vector<double>> m_J_L_q{};
  LinAlg::Matrix<std::vector<double>> m_J_L_q_on_qr{};
  std::vector<double> m_q{};
};

} // namespace SphericalBessel
