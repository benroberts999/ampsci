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
    for (auto L = 0; L <= max_L; ++L) {
      for (auto iq = 0ul; iq < m_J_L_q.cols(); ++iq) {
        const auto tq = q[iq];
        m_J_L_q[std::size_t(L)][iq] = fillBesselVec_kr(L, tq, r);
      }
    }
  }

  //----------------------------------------------------------------------------
  //! Direct access to jL(q*r) for specific q grid index
  const std::vector<double> &at(std::size_t L, std::size_t iq) const {
    return m_J_L_q.atc(std::size_t(L), iq);
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
  std::vector<double> m_q{};
};

} // namespace SphericalBessel
